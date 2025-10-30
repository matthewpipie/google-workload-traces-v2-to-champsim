// Written by Kan Zhu, modified by Matthew Giordano

// DynamoRIO imports
#include "dr_api.h"
#include "drmemtrace/scheduler.h"
// #include "decode_cache.h"
#include "view.h"
#include "instr_api.h"

// ChampSim imports
#include "trace_instruction.h"

// std imports
#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <thread>
#include <vector>
#undef NDEBUG
#include <cassert>
#include <string>
#include <zlib.h>    // For gzopen, gzwrite, gzclose
#include <cstdlib>   // For std::to_string, atoi, rand()
#include <atomic>    // For std::atomic
#include <getopt.h>  // For getopt_long
#include <sys/resource.h>
#include <chrono>
#include <array>
#include <vector>
#include <unordered_map>

// my headers
#include "safe_num_cast.hpp"
#include "threadsafe_rand.hpp"

using namespace dynamorio::drmemtrace;
using namespace champsim;

// -----------------------------------------------------------------------------
// Fault Types and Names (for debugging/statistics)
enum FaultType {
    FAULT_UNKNOWN_INST,
    FAULT_UNKNOWN_CONDITIONAL_BRANCH_TAKEN,
    FAULT_DATA_PC_MISMATCH,
    FAULT_TOO_MANY_SRC_REGISTERS,
    FAULT_TOO_MANY_DST_REGISTERS,
    FAULT_TOO_MANY_SRC_MEMORY,
    FAULT_TOO_MANY_DST_MEMORY,
    FAULT_MAX
};

const char *fault_names[FAULT_MAX + 1] = {
    "Unknown instruction",
    "Unknown conditional branch taken",
    "Data PC mismatch",
    "Too many source registers",
    "Too many destination registers",
    "Too many source memory references",
    "Too many destination memory references",
    "MAX"
};

// Thread-safe fault counters.
std::atomic<uint64_t> failure_counts[FAULT_MAX + 1] = {0};

// -----------------------------------------------------------------------------
// Branch Types and Names (for statistics)
enum BranchType {
    BRANCH_DIRECT_JUMP, // address known
    BRANCH_INDIRECT_JUMP, // address in reg/mem
    BRANCH_DIRECT_CALL, // address known
    BRANCH_INDIRECT_CALL, // address in reg/mem
    BRANCH_RETURN, // address on stack
    BRANCH_TAKEN_JUMP, // only used on conditional taken direct jumps
    BRANCH_UNTAKEN_JUMP, // only used on conditional untaken direct jumps
    BRANCH_CONDITIONAL_JUMP, // depricated, should always be 0 in practice
    BRANCH_INVALID,  // For non-branch or unknown types.
    BRANCH_MAX
};

const char* branch_type_names[] = {
    "Direct Jump",
    "Indirect Jump",
    "Direct Call",
    "Indirect Call",
    "Return",
    "Taken Jump",
    "Untaken Jump",
    "Conditional Jump",
    "Invalid Branch"
};

// Convert from DynamoRIO branch type to our internal branch type
BranchType getBranchType(trace_type_t instr_type) {
    switch (instr_type) {
        case TRACE_TYPE_INSTR_DIRECT_JUMP:      return BRANCH_DIRECT_JUMP;
        case TRACE_TYPE_INSTR_INDIRECT_JUMP:    return BRANCH_INDIRECT_JUMP;
        case TRACE_TYPE_INSTR_DIRECT_CALL:      return BRANCH_DIRECT_CALL;
        case TRACE_TYPE_INSTR_INDIRECT_CALL:    return BRANCH_INDIRECT_CALL;
        case TRACE_TYPE_INSTR_RETURN:           return BRANCH_RETURN;
        case TRACE_TYPE_INSTR_TAKEN_JUMP:       return BRANCH_TAKEN_JUMP;
        case TRACE_TYPE_INSTR_UNTAKEN_JUMP:     return BRANCH_UNTAKEN_JUMP;
        case TRACE_TYPE_INSTR_CONDITIONAL_JUMP: return BRANCH_CONDITIONAL_JUMP;
        default:                                return BRANCH_INVALID;
    }
}

// Global per-branch-type overflow counters.
std::array<std::atomic<uint64_t>, BRANCH_MAX> branch_dst_overflow_counts = {0,0,0,0,0,0,0,0,0};
std::array<std::atomic<uint64_t>, BRANCH_MAX> branch_src_overflow_counts = {0,0,0,0,0,0,0,0,0};

static std::mutex print_mutex;

// -----------------------------------------------------------------------------
// Per-thread statistics structure.
struct SimStats {
    uint64_t total_insts = 0;
    uint64_t branch_count = 0;
    std::array<uint64_t, BRANCH_MAX> branch_type_counts = {0,0,0,0,0,0,0,0,0}; // Indexed by BranchType 0..7.
    ushort lowest_reg_id_seen = 65535;
    ushort highest_reg_id_seen = 0;
};

#define TESTANY(mask, var) (((mask) & (var)) != 0)


// -----------------------------------------------------------------------------
// Helper functions for register assignment.
// For non-branch instructions, we add an offset of 100 to registers;
// for branch instructions, we use the register value as-is.
// Why?
// Unfortunately, the Google traces anonymized all registers and combined several,
// so there is no way to map them to real registers.
// However, ChampSim has 3 "special" registers, where reading/writing to them have special effects on the frontend unlike most general purpose registers:
// the stack pointer, the instruction pointer, and the flags register.
// Namely, for normal instructions, the stack pointer's dependencies are less strict than GPRs, because
// one can almost always determine the new value of the stack pointer post-decode instead of having to wait post-writeback.
// Google anonymized registers, so we don't know which is the stack pointer, so ChampSim won't be able to make this optimization,
// leading to slower traces.
// The other use of the SP and IP in ChampSim is to reverse engineer the branch type.
// Instead of containing branch type in the trace struct,
// ChampSim reverse engineers this via looking at:
// * `is_branch`
// * `branch_taken`
// * `writes_sp`, `reads_sp`
// * `writes_ip`, `reads_ip`
// * `reads_flags` (luckily it seems like the flags register was not anonymized... I hope not, anyways!)
// * `reads_other` (any other register)
// This code is present in ChampSim/inc/instruction.h
// Why is this all relevant?
// ChampSim uses magic register numbers to refer to the SP and IP. 
// DynamoRIO has different magic register numbers for each register.
// We need to make sure we don't accidently create false dependencies with the magic SP and IP register numbers.
// So, we add 100 to each register to completely avoid the SP, IP, and flags registers, which are, in practice, at offset 6, 26, and 25.
// When adding a branch register dependency, we are adding a dependency to SP, IP, or flags, so we don't add anything.
//
// Also, there is one more important caveat.
// ChampSim assumes x86 has a limit of 4 input registers and 4 input memory addresses per instruction,
// and a limit of 2 output registers and 2 output memory address per instruction.
// However, the Google anonymization process seems to, for some reason, occasionally give instructions 5 inputs or 3 outputs.
// In practice, this appears to be super rare (100,000 faults when processing billions of instructions), but it varies per workload,
// so be sure to note it.

constexpr ssize_t REG_DEPENDENCY_OFFSET = 100;

static_assert(REG_DEPENDENCY_OFFSET > REG_STACK_POINTER,
              "Register dependency offset must exceed stack pointer magic register ID");
static_assert(REG_DEPENDENCY_OFFSET > REG_FLAGS,
              "Register dependency offset must exceed flags magic register ID");
static_assert(REG_DEPENDENCY_OFFSET > REG_INSTRUCTION_POINTER,
              "Register dependency offset must exceed instruction pointer magic register ID");

constexpr unsigned char RAND_REG_RANGE_MIN = REG_DEPENDENCY_OFFSET;
constexpr unsigned char RAND_REG_RANGE_MAX = RAND_REG_RANGE_MIN + 15; // 16 registers roughly, this is an inclusive max

constexpr unsigned REG_SRC_NUM_INSTR = NUM_INSTR_SOURCES;
constexpr unsigned REG_DST_NUM_INSTR = NUM_INSTR_DESTINATIONS;
constexpr unsigned MEM_SRC_NUM_INSTR = NUM_INSTR_SOURCES;
constexpr unsigned MEM_DST_NUM_INSTR = NUM_INSTR_DESTINATIONS;

bool addSrcRegister(input_instr &champsim_instr, unsigned &srcCount, ushort reg, bool verbose) {
    if (srcCount >= REG_SRC_NUM_INSTR) {
        if (verbose) std::cerr << "Too many source registers" << std::endl;
        failure_counts[FAULT_TOO_MANY_SRC_REGISTERS]++;
        return false;
    }
    champsim_instr.source_registers[srcCount++] = safe_num_cast<unsigned char>(reg + REG_DEPENDENCY_OFFSET);
    return true;
}

bool addDstRegister(input_instr &champsim_instr, unsigned &dstCount, ushort reg, bool verbose) {
    if (dstCount >= REG_DST_NUM_INSTR) {
        if (verbose) std::cerr << "Too many destination registers" << std::endl;
        failure_counts[FAULT_TOO_MANY_DST_REGISTERS]++;
        return false;
    }
    champsim_instr.destination_registers[dstCount++] = safe_num_cast<unsigned char>(reg + REG_DEPENDENCY_OFFSET);
    return true;
}

// Branch register helpers: no offset added.
bool addSrcRegisterForBranch(input_instr &champsim_instr, unsigned &srcCount, ushort reg, bool verbose, BranchType bType) {
    if (srcCount >= REG_SRC_NUM_INSTR) {
        if (verbose) std::cerr << "Too many source registers for branch type " 
                              << branch_type_names[bType] << std::endl;
        branch_src_overflow_counts[bType]++;
        return false;
    }
    champsim_instr.source_registers[srcCount++] = safe_num_cast<unsigned char>(reg);
    return true;
}

bool addDstRegisterForBranch(input_instr &champsim_instr, unsigned &dstCount, ushort reg, bool verbose, BranchType bType) {
    if (dstCount >= REG_DST_NUM_INSTR) {
        if (verbose) std::cerr << "Too many destination registers for branch type " 
                              << branch_type_names[bType] << std::endl;
        branch_dst_overflow_counts[bType]++;
        return false;
    }
    champsim_instr.destination_registers[dstCount++] = safe_num_cast<unsigned char>(reg);
    return true;
}

bool addSrcRegisters(input_instr &champsim_instr, unsigned &srcCount, const std::vector<ushort> &regs, bool verbose) {
    for (ushort reg : regs) {
        if (!addSrcRegister(champsim_instr, srcCount, reg, verbose))
            return false;
    }
    return true;
}

bool addDstRegisters(input_instr &champsim_instr, unsigned &dstCount, const std::vector<ushort> &regs, bool verbose) {
    for (ushort reg : regs) {
        if (!addDstRegister(champsim_instr, dstCount, reg, verbose))
            return false;
    }
    return true;
}

bool addSrcRegistersForBranch(input_instr &champsim_instr, unsigned &srcCount, const std::vector<ushort> &regs, bool verbose, BranchType bType) {
    for (ushort reg : regs) {
        if (!addSrcRegisterForBranch(champsim_instr, srcCount, reg, verbose, bType))
            return false;
    }
    return true;
}

bool addDstRegistersForBranch(input_instr &champsim_instr, unsigned &dstCount, const std::vector<ushort> &regs, bool verbose, BranchType bType) {
    for (ushort reg : regs) {
        if (!addDstRegisterForBranch(champsim_instr, dstCount, reg, verbose, bType))
            return false;
    }
    return true;
}

// -----------------------------------------------------------------------------
// Helper: Consolidated branch-specific register assignments.
// Uses specialized branch helpers (without offset).
bool assignBranchRegisters(input_instr &champsim_input_instr, unsigned &srcCount, unsigned &dstCount, trace_type_t branchType, bool verbose) {
    BranchType bType = getBranchType(branchType);
    switch (branchType) {
        case TRACE_TYPE_INSTR_DIRECT_JUMP:
            // No extra register needed.
            // Caller already adds IP as a dest register.
            break;
        case TRACE_TYPE_INSTR_INDIRECT_JUMP:
            // An indirect branch should really depend on a register
            // If it doesn't, we should artificially add one to make sure ChampSim recognizes it as an indirect branch
            // This could be considered a fault, I suppose
            if (srcCount == 0) {
                if (verbose) std::cout << "Adding random register for indirect jump" << std::endl;
                if (!addSrcRegisterForBranch(champsim_input_instr, srcCount, threadsafe_rand<ushort>(RAND_REG_RANGE_MIN, RAND_REG_RANGE_MAX), verbose, bType))
                    return false;
            }
            break;
        case TRACE_TYPE_INSTR_CONDITIONAL_JUMP:
        case TRACE_TYPE_INSTR_TAKEN_JUMP:
        case TRACE_TYPE_INSTR_UNTAKEN_JUMP:
            // Conditional branches read and write IP, and read (flags || other)
            // Similarly to above, generate a random dependency if none exists
            if (srcCount == 0) {
                if (verbose) std::cout << "Adding random register for conditional jump" << std::endl;
                if (threadsafe_rand<int>(0, 1) == 0) {
                    // random register
                    if (!addSrcRegisterForBranch(champsim_input_instr, srcCount, threadsafe_rand<ushort>(RAND_REG_RANGE_MIN, RAND_REG_RANGE_MAX), verbose, bType))
                        return false;
                } else {
                    // flags
                    if (!addSrcRegisterForBranch(champsim_input_instr, srcCount, REG_FLAGS, verbose, bType))
                        return false;
                }
            }
            if (!addSrcRegisterForBranch(champsim_input_instr, srcCount, REG_INSTRUCTION_POINTER, verbose, bType))
                return false;
            break;
        case TRACE_TYPE_INSTR_DIRECT_CALL:
            // Reads IP, SP, writes IP, SP, doesn't read flags or other
            // To do this, we overwrite all src and dest regs
            srcCount = 0;
            if (!addSrcRegistersForBranch(champsim_input_instr, srcCount, {REG_STACK_POINTER, REG_INSTRUCTION_POINTER}, verbose, bType))
                return false;
            dstCount = 0;
            if (!addDstRegistersForBranch(champsim_input_instr, dstCount, {REG_STACK_POINTER, REG_INSTRUCTION_POINTER}, verbose, bType))
                return false;
            // Zero out remaining fields
            for (unsigned i = srcCount; i < REG_SRC_NUM_INSTR; i++) {
                if (!addSrcRegisterForBranch(champsim_input_instr, srcCount, 0, verbose, bType)) return false;
            }
            for (unsigned i = dstCount; i < REG_DST_NUM_INSTR; i++) {
                if (!addDstRegisterForBranch(champsim_input_instr, dstCount, 0, verbose, bType)) return false;
            }
            break;
        case TRACE_TYPE_INSTR_INDIRECT_CALL:
            // Reads IP, SP
            // Writes IP, SP
            // Reads other, doesn't read flags
            if (srcCount == 0) {
                if (verbose) std::cout << "Adding random register for indirect call" << std::endl;
                if (!addSrcRegisterForBranch(champsim_input_instr, srcCount, threadsafe_rand<ushort>(RAND_REG_RANGE_MIN, RAND_REG_RANGE_MAX), verbose, bType))
                    return false;
            }
            if (!addSrcRegistersForBranch(champsim_input_instr, srcCount, {REG_STACK_POINTER, REG_INSTRUCTION_POINTER}, verbose, bType))
                return false;
            if (!addDstRegisterForBranch(champsim_input_instr, dstCount, REG_STACK_POINTER, verbose, bType))
                return false;
            break;
        case TRACE_TYPE_INSTR_RETURN:
            // Reads SP, not IP
            // writes SP, IP
            if (!addSrcRegisterForBranch(champsim_input_instr, srcCount, REG_STACK_POINTER, verbose, bType))
                return false;
            if (!addDstRegisterForBranch(champsim_input_instr, dstCount, REG_STACK_POINTER, verbose, bType))
                return false;
            break;
        default:
            break;
    }
    return true;
}

void update_reg_bounds(ushort reg, SimStats &stats) {
    if (reg < stats.lowest_reg_id_seen) {
        stats.lowest_reg_id_seen = reg;
    }
    if (reg > stats.highest_reg_id_seen) {
        stats.highest_reg_id_seen = reg;
    }
}

// -----------------------------------------------------------------------------
// update_inst_registers: Update input_instr registers for the given instruction.
void update_inst_registers(void *dcontext, memref_t record, instr_t dr_instr, input_instr &champsim_input_instr, SimStats &stats, bool verbose = false) {
    (void)(dcontext); // "use" it
    unsigned srcCount = 0;
    unsigned dstCount = 0;
    uint used_flag = instr_get_arith_flags(&dr_instr, DR_QUERY_DEFAULT);

    for (uint i = 0; i < safe_num_cast<uint>(instr_num_srcs(&dr_instr)); i++) {
        opnd_t opnd = instr_get_src(&dr_instr, i);
        for (int opnum = 0; opnum < opnd_num_regs_used(opnd); opnum++) {
            reg_id_t reg = opnd_get_reg_used(opnd, opnum);
            update_reg_bounds(reg, stats);
            if (verbose) {
                std::cout << "src register " << i << "," << opnum << ": " << reg << std::endl;
            }
            if (!addSrcRegister(champsim_input_instr, srcCount, safe_num_cast<ushort>(reg), verbose))
                return;
        }
    }

    for (uint i = 0; i < safe_num_cast<uint>(instr_num_dsts(&dr_instr)); i++) {
        opnd_t opnd = instr_get_dst(&dr_instr, i);
        for (int opnum = 0; opnum < opnd_num_regs_used(opnd); opnum++) {
            reg_id_t reg = opnd_get_reg_used(opnd, opnum);
            update_reg_bounds(reg, stats);
            if (verbose) {
                std::cout << "dst register " << i << "," << opnum << ": " << reg << std::endl;
            }
            if (!addDstRegister(champsim_input_instr, dstCount, safe_num_cast<ushort>(reg), verbose))
                return;
        }
    }

    if (TESTANY(EFLAGS_WRITE_ARITH, used_flag)) {
        if (verbose) std::cout << "EFLAGS is written" << std::endl;
        if (!addDstRegisterForBranch(champsim_input_instr, dstCount, REG_FLAGS, verbose, BRANCH_INVALID))
            return;
    }
    if (TESTANY(EFLAGS_READ_ARITH, used_flag)) {
        if (verbose) std::cout << "EFLAGS is read" << std::endl;
        if (!addSrcRegisterForBranch(champsim_input_instr, srcCount, REG_FLAGS, verbose, BRANCH_INVALID))
            return;
    }

    if (record.instr.type == TRACE_TYPE_INSTR) {
        return;
    } else {
        BranchType bType = getBranchType(record.instr.type);
        dstCount = 0;
        // IP is always written on a branch
        if (!addDstRegisterForBranch(champsim_input_instr, dstCount, REG_INSTRUCTION_POINTER, verbose, bType)) {
            return;
        }
    }

    if (!assignBranchRegisters(champsim_input_instr, srcCount, dstCount, record.instr.type, verbose)) {
        return;
    }
}

// -----------------------------------------------------------------------------
// update_branch_info: Update branch flags and fault counts.
void update_branch_info(memref_t record, input_instr &champsim_input_instr) {
    switch (record.instr.type) {
        case TRACE_TYPE_INSTR:
            champsim_input_instr.is_branch = false;
            champsim_input_instr.branch_taken = false;
            break;
        case TRACE_TYPE_INSTR_DIRECT_JUMP:
        case TRACE_TYPE_INSTR_INDIRECT_JUMP:
        case TRACE_TYPE_INSTR_DIRECT_CALL:
        case TRACE_TYPE_INSTR_INDIRECT_CALL:
        case TRACE_TYPE_INSTR_RETURN:
        case TRACE_TYPE_INSTR_TAKEN_JUMP:
            champsim_input_instr.is_branch = true;
            champsim_input_instr.branch_taken = true;
            break;
        case TRACE_TYPE_INSTR_UNTAKEN_JUMP:
            champsim_input_instr.is_branch = true;
            champsim_input_instr.branch_taken = false;
            break;
        case TRACE_TYPE_INSTR_CONDITIONAL_JUMP:
            champsim_input_instr.is_branch = true;
            champsim_input_instr.branch_taken = false;
            std::cerr << "unknown conditional jump" << std::endl;
            failure_counts[FAULT_UNKNOWN_CONDITIONAL_BRANCH_TAKEN]++;
            break;
        default:
            champsim_input_instr.is_branch = false;
            champsim_input_instr.branch_taken = false;
            std::cerr << "Unknown instruction type: " << record.instr.type << std::endl;
            failure_counts[FAULT_UNKNOWN_INST]++;
            break;
    }
}

struct thread_args {
    void* dcontext; // DynamoRIO context for decoding, note: not fully thread safe, but should be "safe enough" for our usage
    scheduler_t::stream_t *stream; // Instruction stream
    const int thread_id;
    const bool verbose;
    SimStats &stats;
    const std::string &output_file_without_suffix;
    const uint64_t instruction_cap_global; // If nonzero, stop as soon as all threads combined have processed this many instructions
    std::atomic<size_t> &g_global_inst_count; // number of instructions processed, across all threads, if above is nonzero
    const uint64_t instruction_cap_local; // If nonzero, stop this thread once this number of instructions have been reached
};

// -----------------------------------------------------------------------------
// simulate_core: Processes records from a given stream.
void simulate_core(thread_args &args) {

    std::unordered_map<size_t, instr_t> pc_instr_map;

    auto start_time = std::chrono::high_resolution_clock::now();

    int output_file_zero_width = 4; // number of digits you want, e.g. 0001, 0002, ...
    std::ostringstream oss;
    oss << std::setw(output_file_zero_width) << std::setfill('0') << args.thread_id;

    std::string filename = args.output_file_without_suffix + "_" + oss.str() + ".champsim.gz";
    
    gzFile gz_out = gzopen(filename.c_str(), "w");
    if (!gz_out) {
        std::cerr << "Failed to open gz file " << filename << "\n";
        return;
    }

    // aliased because we use it so much
    scheduler_t::stream_t *stream = args.stream;
    SimStats &stats = args.stats;
    
    if (args.verbose) {
        std::cout << "Thread " << args.thread_id << " processing filetype " 
                  << stream->get_filetype() << "\n";
    }
    
    memref_t record;
    scheduler_t::stream_status_t status = stream->next_record(record);
    
    while (status != scheduler_t::STATUS_EOF) {
        if (status == scheduler_t::STATUS_WAIT || status == scheduler_t::STATUS_IDLE) {
            std::this_thread::yield();
            status = stream->next_record(record);
            continue;
        }
        assert(status == scheduler_t::STATUS_OK);
        
        if (type_is_instr(record.instr.type)) {
            stats.total_insts++;
            if (stats.total_insts % 1000000 == 0) {
                auto elapsed_time = std::chrono::high_resolution_clock::now() - start_time;
                double elapsed_seconds = std::chrono::duration<double>(elapsed_time).count();
                print_mutex.lock();
                std::cout << " Core " <<  args.thread_id << " Processed " << stats.total_insts << " instructions, "
                          << " Elapsed time: "
                          << elapsed_seconds << " seconds" << ", MI/s: "
                          << ((double)(stats.total_insts) / 1000000.0) / elapsed_seconds << "\n";
                print_mutex.unlock();
            }

            if (args.verbose) {
                std::cout << std::left << std::setw(12) << "ifetch" << std::right
                          << std::setw(2) << record.instr.size << " byte(s) @ 0x" << std::hex
                          << std::setfill('0') << std::setw(sizeof(void *) * 2) << record.instr.addr
                          << std::dec << std::setfill(' ') << "\n";
            }

            // We need to "decode" the instruction to get all of its operands, etc.
            // To make this faster, we cache known decodes in the map `pc_instr_map`
            size_t pc = record.instr.addr;
            auto it = pc_instr_map.find(pc);
            instr_t dr_instr;
            if (it != pc_instr_map.end() && !record.instr.encoding_is_new) {
                // cache hit
                dr_instr = it->second;
            } else {
                // cache miss
                // decode
                instr_init(args.dcontext, &dr_instr);
                const app_pc decode_pc = reinterpret_cast<app_pc>(pc);
                // It's unclear if we need to use decode_from_copy() or if we can just use decode()
                // It's also unclear if the second argument here is correct
                app_pc _nextpc = decode_from_copy(args.dcontext, record.instr.encoding, decode_pc, &dr_instr);
                (void)(_nextpc);
                // Write into decode cache
                pc_instr_map[pc] = dr_instr;
            }

            // Now, convert DynamoRIO instruction into ChampSim instruction
            input_instr champsim_inst;
            memset(&champsim_inst, 0, sizeof(champsim_inst)); // zero out all fields
            // There are four parts to a CS instruction:
            //  * ip
            //  * branch info
            //  * src, dest regs
            //  * src, dest mems
            
            // Copy src, dest regs
            update_inst_registers(args.dcontext, record, dr_instr, champsim_inst, stats, args.verbose);
            
            // Copy ip
            champsim_inst.ip = record.instr.addr;
            assert(champsim_inst.ip != 0);

            // Copy branch info
            update_branch_info(record, champsim_inst);
            if (champsim_inst.is_branch) {
                stats.branch_count++;
                BranchType bType = getBranchType(record.instr.type);
                stats.branch_type_counts[bType]++;
            }
            
            // Copy src, dest memory dependencies
            // In order to get memory source and destination information about this instruction, 
            // we need to read the next record.
            // Sort of annoyingly repetitive but alas...
            memref_t new_record;
            status = stream->next_record(new_record);
            size_t inst_source_memory = 0;
            size_t inst_dest_memory = 0;
            while (status != scheduler_t::STATUS_EOF &&
                   (new_record.instr.type == TRACE_TYPE_READ || new_record.instr.type == TRACE_TYPE_WRITE)) {
                if (status == scheduler_t::STATUS_WAIT || status == scheduler_t::STATUS_IDLE) {
                    std::this_thread::yield();
                    status = stream->next_record(new_record);
                    continue;
                }
                assert(status == scheduler_t::STATUS_OK);
                if (new_record.data.pc != record.instr.addr) {
                    std::cerr << "Memory access does not match instruction address\n";
                    failure_counts[FAULT_DATA_PC_MISMATCH]++;
                }
                if (new_record.instr.type == TRACE_TYPE_READ) {
                    if (inst_source_memory >= MEM_SRC_NUM_INSTR) {
                        if (args.verbose) std::cerr << "Too many source memory" << std::endl;
                        failure_counts[FAULT_TOO_MANY_SRC_MEMORY]++;
                    } else {
                        champsim_inst.source_memory[inst_source_memory++] = new_record.instr.addr;
                    }
                    assert(new_record.instr.addr != 0);
                } else {
                    assert(new_record.instr.type == TRACE_TYPE_WRITE);
                    if (inst_dest_memory >= MEM_DST_NUM_INSTR) {
                        if (args.verbose) std::cerr << "Too many destination memory" << std::endl;
                        failure_counts[FAULT_TOO_MANY_DST_MEMORY]++;
                    } else {
                        champsim_inst.destination_memory[inst_dest_memory++] = new_record.instr.addr;
                    }
                    assert(new_record.instr.addr != 0);
                }
                status = stream->next_record(new_record);
            }

            // Done, write out to gz file
            gzwrite(gz_out, &champsim_inst, sizeof(champsim_inst));
            record = new_record; // copy over the next record from the above code block

            // Break conditions
            if (args.instruction_cap_local != 0 && stats.total_insts >= args.instruction_cap_local) {
                break;
            }
            if (args.instruction_cap_global != 0) {
                uint64_t new_val = args.g_global_inst_count.fetch_add(1) + 1;
                if (new_val >= args.instruction_cap_global) {
                    break;
                }
            }

        } else {
            status = stream->next_record(record);
        }
    }
    gzclose(gz_out);
}

// -----------------------------------------------------------------------------
// run_scheduler: Spawns one thread per core and writes output files.
void run_scheduler(const std::string &trace_directory, bool verbose,
                   std::vector<SimStats> &all_stats,
                   const std::string &output_file_path,
                   const std::string &output_file_name,
                   int num_cores,
                   uint64_t instruction_cap_global,
                   uint64_t instruction_cap_local
                ) {
    scheduler_t scheduler;
    std::vector<scheduler_t::input_workload_t> sched_inputs;
    sched_inputs.emplace_back(trace_directory);
    
    scheduler_t::scheduler_options_t sched_ops(scheduler_t::MAP_TO_ANY_OUTPUT,
                                               scheduler_t::DEPENDENCY_TIMESTAMPS,
                                               scheduler_t::SCHEDULER_DEFAULTS);
    scheduler_t::scheduler_status_t status = scheduler.init(sched_inputs, num_cores, std::move(sched_ops));
    if (status != scheduler_t::STATUS_SUCCESS) {
            std::cout << "error: " << status << std::endl;
            std::cout << "estring: " << scheduler.get_error_string() << std::endl;
        assert(false);
    }
    
    void *dcontext = dr_standalone_init();
    auto filetype_ = scheduler.get_stream(0)->get_filetype();
    
    if (TESTANY(OFFLINE_FILE_TYPE_ARCH_ALL & ~OFFLINE_FILE_TYPE_ARCH_REGDEPS, filetype_) &&
        !TESTANY(build_target_arch_type(), safe_num_cast<int>(filetype_))) {
        std::string error_string_ = std::string("Architecture mismatch: trace recorded on ") +
                                    trace_arch_string(static_cast<offline_file_type_t>(filetype_)) +
                                    " but tool built for " +
                                    trace_arch_string(build_target_arch_type());
        std::cerr << error_string_ << std::endl;
    }
    
    if (TESTANY(OFFLINE_FILE_TYPE_ARCH_REGDEPS, filetype_)) {
        dr_set_isa_mode(dcontext, DR_ISA_REGDEPS, nullptr);
    }
    
    dr_disasm_flags_t flags = IF_X86_ELSE(
        DR_DISASM_ATT,
        IF_AARCH64_ELSE(DR_DISASM_DR, IF_RISCV64_ELSE(DR_DISASM_RISCV, DR_DISASM_ARM)));
    if (TESTANY(OFFLINE_FILE_TYPE_ARCH_REGDEPS, filetype_)) {
        flags = DR_DISASM_DR;
    }
    disassemble_set_syntax(flags);

    std::string output_better = output_file_path;
    if (!output_better.empty() && output_better.back() != '/' && output_better.back() != '\\') {
        output_better += "/";
    }
    output_better += output_file_name;

    all_stats.resize(safe_num_cast<unsigned>(num_cores));
    std::vector<std::thread> threads;
    std::vector<thread_args> args;
    threads.reserve(safe_num_cast<unsigned>(num_cores));
    args.reserve(safe_num_cast<unsigned>(num_cores));
    std::atomic<uint64_t> g_global_inst_count{0};
    for (int i = 0; i < num_cores; ++i) {
        args.push_back(thread_args {
            /* .dcontext = */ dcontext,
            /* .stream = */ scheduler.get_stream(i),
            /* .thread_id = */ i,
            /* .verbose = */ verbose,
            /* .stats = */ all_stats[safe_num_cast<unsigned>(i)],
            /* .output_file_without_suffix = */ output_better,
            /* .instruction_cap_global = */ instruction_cap_global,
            /* .g_global_inst_count = */ std::ref(g_global_inst_count),
            /* .instruction_cap_local = */ instruction_cap_local,
        });
    }
    // Must do this in a separate for loop to avoid race condition of vector getting resized and moving elements before thread starts
    for (size_t i = 0; i < safe_num_cast<size_t>(num_cores); ++i) {
        thread_args &my_args = args[i];
        threads.emplace_back([&my_args]() {
            simulate_core(my_args);
        });
    }
    for (std::thread &thread : threads) {
        thread.join();
    }
}


void setrlim() {
    rlimit lim;
    
    if (getrlimit(RLIMIT_NOFILE, &lim) != 0) {
        perror("getrlimit");
        exit(1);
    }

    std::cout << "Current limits:\n";
    std::cout << "  Soft: " << (lim.rlim_cur == RLIM_INFINITY ? "infinity" : std::to_string(lim.rlim_cur)) << "\n";
    std::cout << "  Hard: " << (lim.rlim_max == RLIM_INFINITY ? "infinity" : std::to_string(lim.rlim_max)) << "\n";

    lim.rlim_cur = lim.rlim_max;

    if (setrlimit(RLIMIT_NOFILE, &lim) != 0) {
        perror("setrlimit");
        exit(1);
    }

    std::cout << "RLIMIT_NOFILE set to " << lim.rlim_cur << ".\n";
}

// -----------------------------------------------------------------------------
// Main: Parse command-line arguments, run the scheduler, and print statistics.
int main(int argc, char *argv[]) {
    std::string trace_directory;
    std::string output_file_path = ".";
    std::string output_file_name = "bravo";
    bool verbose = false;
    uint64_t percore_target_count = 0;
    uint64_t global_target_count = 0;
    int num_cores = 0;

    int opt;
    int option_index = 0;
    static struct option long_options[] = {
        {"trace_folder",     required_argument, 0, 't'},
        {"output_file_path", required_argument, 0, 'p'},
        {"output_file_name", required_argument, 0, 'f'},  // Changed from -n to -f.
        {"num_cores",        required_argument, 0, 'n'},  // New option for number of cores.
        {"verbose",            no_argument,       0, 'v'},
        {"percore_target_count",    required_argument, 0, 'c'}, // cap the per-core instruction count
        {"global_target_count",     required_argument, 0, 'g'}, // cap the global instruction count
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "t:p:f:n:qc:g:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 't':
                trace_directory = optarg;
                break;
            case 'p':
                output_file_path = optarg;
                break;
            case 'f':  // File name option.
                output_file_name = optarg;
                break;
            case 'n':  // Number of cores option.
                num_cores = std::stoi(optarg);
                break;
            case 'v':
                verbose = true;
                break;
            case 'c':
                percore_target_count = std::stoull(optarg);
                break;
            case 'g':
                global_target_count = std::stoull(optarg);
                break;
            default:
                std::cerr << "Usage: " << argv[0]
                          << " --trace_folder <folder> [--output_file_path <path>] [--output_file_name <name>] [--num_cores <num>] [--verbose] [--percore_target_count <num>] [--global_target_count <num>]\n";
                return 1;
        }
    }

    if (trace_directory.empty()) {
        std::cerr << "Error: --trace_folder is required.\n";
        std::cerr << "Usage: " << argv[0]
                  << " --trace_folder <folder> [--output_file_path <path>] [--output_file_name <name>] [--num_cores <num>] [--verbose] [--percore_target_count <num>] [--global_target_count <num>]\n";
        return 1;
    }

    setrlim();

    std::vector<SimStats> thread_stats;
    auto start = std::chrono::high_resolution_clock::now();
    run_scheduler(trace_directory, verbose, thread_stats, output_file_path, output_file_name, num_cores, global_target_count, percore_target_count);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    uint64_t total_insts = 0;
    uint64_t total_branches = 0;
    std::array<uint64_t, 8> total_branch_types = {0, 0, 0, 0, 0, 0, 0, 0};
    for (const auto &s : thread_stats) {
        total_insts += s.total_insts;
        total_branches += s.branch_count;
        for (size_t i = 0; i < 8; i++) {
            total_branch_types[i] += s.branch_type_counts[i];
        }
    }
    
    std::cout << "\n=== Simulation Statistics ===\n";
    std::cout << "Total instructions processed: " << total_insts << "\n";
    std::cout << "Total branch instructions:    " << total_branches << "\n";
    
    std::cout << "\nBranch Types for each core:\n";
    for (size_t i = 0; i < thread_stats.size(); i++) {
        std::cout << "Thread " << i << ":\n";
        for (size_t j = 0; j < 8; j++) {
            std::cout << branch_type_names[j] << ": " << thread_stats[i].branch_type_counts[j] << "\n";
        }
    }

    std::cout << "\nFaults:\n";
    for (size_t i = 0; i < FAULT_MAX; ++i) {
        std::cout << fault_names[i] << ": " << failure_counts[i] << "\n";
    }
    
    std::cout << "\nBranch Destination Register Overflows:\n";
    for (size_t i = 0; i < 8; i++) {
        std::cout << branch_type_names[i] << ": " << branch_dst_overflow_counts[i] << "\n";
    }
    
    std::cout << "\nBranch Source Register Overflows:\n";
    for (size_t i = 0; i < 8; i++) {
        std::cout << branch_type_names[i] << ": " << branch_src_overflow_counts[i] << "\n";
    }
    
    for (size_t i = 0; i < thread_stats.size(); i++) {
        std::cout << "Thread " << i << " processed " << thread_stats[i].total_insts
        << " instructions, reg bounds = (" << thread_stats[i].lowest_reg_id_seen << "," << thread_stats[i].highest_reg_id_seen << ")\n";
    }
    
    std::cout << "\nTime taken: " << elapsed.count() << " seconds\n";
    std::cout << "Seconds per million instructions: " << elapsed.count() / ((double)(total_insts) / 1e6) << "\n";

    return 0;
}
