// Written originally by Kan Zhu
// Modified by Matthew Giordano

#include <dr_api.h>
#include "drmemtrace/scheduler.h"
#include "decode_cache.h"
#include "view.h"

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
#include <trace_instruction.h>
#include <zlib.h>    // For gzopen, gzwrite, gzclose
#include <cstdlib>   // For std::to_string, atoi, rand()
#include <atomic>    // For std::atomic
#include <getopt.h>  // For getopt_long
#include "instr_api.h"
#include <sys/resource.h>

#include <chrono>
#include <array>
#include <vector>

#include <unordered_map>

using namespace dynamorio::drmemtrace;
using namespace champsim;

// -----------------------------------------------------------------------------
// Global atomic counters for instruction processing.
std::atomic<size_t> g_global_inst_count(0); // number of instructions processed, across all threads
std::atomic<size_t> g_target_inst_count(0); // number of instructions to stop at, not changed once threads start

// -----------------------------------------------------------------------------
// Fault Types and Names (for debugging/statistics)
enum FaultType {
    FAULT_UNKNOWN_INST,
    FAULT_UNKNOWN_CONDITIONAL_BRANCH_TAKEN,
    FAULT_DATA_PC_MISMATCH,
    FAULT_TOO_MANY_SRC_REGISTERS,
    FAULT_TOO_MANY_DST_REGISTERS,
    FAULT_MAX
};

const char *fault_names[FAULT_MAX + 1] = {
    "Unknown instruction",
    "Unknown conditional branch taken",
    "Data PC mismatch",
    "Too many source registers",
    "Too many destination registers",
    "MAX"
};

// Thread-safe fault counters.
std::atomic<int> failure_counts[FAULT_MAX + 1] = {0};

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
    BRANCH_INVALID  // For non-branch or unknown types.
};

const char* branch_type_names[] = {
    "Direct Jump",
    "Indirect Jump",
    "Direct Call",
    "Indirect Call",
    "Return",
    "Taken Jump",
    "Untaken Jump",
    "Conditional Jump"
};

// Convert from DynamoRIO branch type to our internal branch type
BranchType getBranchType(int instr_type) {
    switch (instr_type) {
        case TRACE_TYPE_INSTR_DIRECT_JUMP:    return BRANCH_DIRECT_JUMP;
        case TRACE_TYPE_INSTR_INDIRECT_JUMP:    return BRANCH_INDIRECT_JUMP;
        case TRACE_TYPE_INSTR_DIRECT_CALL:      return BRANCH_DIRECT_CALL;
        case TRACE_TYPE_INSTR_INDIRECT_CALL:    return BRANCH_INDIRECT_CALL;
        case TRACE_TYPE_INSTR_RETURN:           return BRANCH_RETURN;
        case TRACE_TYPE_INSTR_TAKEN_JUMP:       return BRANCH_TAKEN_JUMP;
        case TRACE_TYPE_INSTR_UNTAKEN_JUMP:     return BRANCH_UNTAKEN_JUMP;
        case TRACE_TYPE_INSTR_CONDITIONAL_JUMP: return BRANCH_CONDITIONAL_JUMP;
        default:                              return BRANCH_INVALID;
    }
}

// Global per-branch-type overflow counters.
std::array<std::atomic<int>, 8> branch_dst_overflow_counts = {0,0,0,0,0,0,0,0};
std::array<std::atomic<int>, 8> branch_src_overflow_counts = {0,0,0,0,0,0,0,0};

static std::mutex print_mutex;

// -----------------------------------------------------------------------------
// Per-thread statistics structure.
struct SimStats {
    size_t total_insts = 0;
    size_t branch_count = 0;
    std::array<size_t, 8> branch_type_counts = {0,0,0,0,0,0,0,0}; // Indexed by BranchType 0..7.
};

// -----------------------------------------------------------------------------
// disasm_info_t: Holds disassembly information.
class disasm_info_t : public decode_info_base_t {
public:
    std::string disasm_;
    // Debug fields could be enabled if desired.
    instr_t *instr = nullptr;
private:
    std::string set_decode_info_derived(void *dcontext, const _memref_instr_t &memref_instr,
                                          instr_t *instr, app_pc decode_pc) {
        this->instr = instr;
        return "";
    }
};

#define TESTANY(mask, var) (((mask) & (var)) != 0)

// -----------------------------------------------------------------------------
// Helper: Initialize decode_cache and set disassembly flags.
std::unique_ptr<decode_cache_t<disasm_info_t>>
init_decode_cache_and_set_flags(void* dcontext, scheduler_t::stream_t* stream) {
    auto decode_cache_ = std::make_unique<decode_cache_t<disasm_info_t>>(
        dcontext, true, true);
    auto filetype_ = stream->get_filetype();
    std::string error_string_ = decode_cache_->init(static_cast<offline_file_type_t>(filetype_));
    return decode_cache_;
}

// -----------------------------------------------------------------------------
// Helper functions for register assignment.
// For non-branch instructions, we add an offset of 100 to registers;
// for branch instructions, we use the register value as-is.

bool addSrcRegister(input_instr &instr, int &srcCount, uint8_t reg, bool verbose) {
    if (srcCount >= 4) {
        if (verbose) std::cerr << "Too many source registers" << std::endl;
        failure_counts[FAULT_TOO_MANY_SRC_REGISTERS]++;
        return false;
    }
    instr.source_registers[srcCount++] = reg + 100;
    return true;
}

bool addDstRegister(input_instr &instr, int &dstCount, uint8_t reg, bool verbose) {
    if (dstCount >= 2) {
        if (verbose) std::cerr << "Too many destination registers" << std::endl;
        failure_counts[FAULT_TOO_MANY_DST_REGISTERS]++;
        return false;
    }
    instr.destination_registers[dstCount++] = reg + 100;
    return true;
}

// Branch register helpers: no offset added.
bool addSrcRegisterForBranch(input_instr &instr, int &srcCount, uint8_t reg, bool verbose, BranchType bType) {
    if (srcCount >= 4) {
        if (verbose) std::cerr << "Too many source registers for branch type " 
                              << branch_type_names[bType] << std::endl;
        branch_src_overflow_counts[bType]++;
        return false;
    }
    instr.source_registers[srcCount++] = reg;
    return true;
}

bool addDstRegisterForBranch(input_instr &instr, int &dstCount, uint8_t reg, bool verbose, BranchType bType) {
    if (dstCount >= 2) {
        if (verbose) std::cerr << "Too many destination registers for branch type " 
                              << branch_type_names[bType] << std::endl;
        branch_dst_overflow_counts[bType]++;
        return false;
    }
    instr.destination_registers[dstCount++] = reg;
    return true;
}

bool addSrcRegisters(input_instr &instr, int &srcCount, const std::vector<uint8_t> &regs, bool verbose) {
    for (uint8_t reg : regs) {
        if (!addSrcRegister(instr, srcCount, reg, verbose))
            return false;
    }
    return true;
}

bool addDstRegisters(input_instr &instr, int &dstCount, const std::vector<uint8_t> &regs, bool verbose) {
    for (uint8_t reg : regs) {
        if (!addDstRegister(instr, dstCount, reg, verbose))
            return false;
    }
    return true;
}

bool addSrcRegistersForBranch(input_instr &instr, int &srcCount, const std::vector<uint8_t> &regs, bool verbose, BranchType bType) {
    for (uint8_t reg : regs) {
        if (!addSrcRegisterForBranch(instr, srcCount, reg, verbose, bType))
            return false;
    }
    return true;
}

bool addDstRegistersForBranch(input_instr &instr, int &dstCount, const std::vector<uint8_t> &regs, bool verbose, BranchType bType) {
    for (uint8_t reg : regs) {
        if (!addDstRegisterForBranch(instr, dstCount, reg, verbose, bType))
            return false;
    }
    return true;
}

// -----------------------------------------------------------------------------
// Helper: Consolidated branch-specific register assignments.
// Uses specialized branch helpers (without offset).
bool assignBranchRegisters(input_instr &input_instr, int &srcCount, int &dstCount, int branchType, bool verbose) {
    BranchType bType = getBranchType(branchType);
    switch (branchType) {
        case TRACE_TYPE_INSTR_DIRECT_JUMP:
            // No extra register needed.
            break;
        case TRACE_TYPE_INSTR_INDIRECT_JUMP:
            if (srcCount == 0) {
                if (verbose) std::cout << "Adding random register for indirect jump" << std::endl;
                if (!addSrcRegisterForBranch(input_instr, srcCount, rand() % 16 + 50, verbose, bType))
                    return false;
            }
            break;
        case TRACE_TYPE_INSTR_CONDITIONAL_JUMP:
        case TRACE_TYPE_INSTR_TAKEN_JUMP:
        case TRACE_TYPE_INSTR_UNTAKEN_JUMP:
            if (!addSrcRegisterForBranch(input_instr, srcCount, REG_INSTRUCTION_POINTER, verbose, bType))
                return false;
            break;
        case TRACE_TYPE_INSTR_DIRECT_CALL:
            srcCount = 0;
            if (!addSrcRegistersForBranch(input_instr, srcCount, {REG_STACK_POINTER, REG_INSTRUCTION_POINTER}, verbose, bType))
                return false;
            dstCount = 0;
            if (!addDstRegistersForBranch(input_instr, dstCount, {REG_STACK_POINTER, REG_INSTRUCTION_POINTER}, verbose, bType))
                return false;
            break;
        case TRACE_TYPE_INSTR_INDIRECT_CALL:
            if (srcCount == 0) {
                if (verbose) std::cout << "Adding random register for indirect call" << std::endl;
                if (!addSrcRegisterForBranch(input_instr, srcCount, rand() % 16 + 50, verbose, bType))
                    return false;
            }
            if (!addSrcRegistersForBranch(input_instr, srcCount, {REG_STACK_POINTER, REG_INSTRUCTION_POINTER}, verbose, bType))
                return false;
            dstCount = 0;
            if (!addDstRegistersForBranch(input_instr, dstCount, {REG_STACK_POINTER, REG_INSTRUCTION_POINTER}, verbose, bType))
                return false;
            break;
        case TRACE_TYPE_INSTR_RETURN:
            if (!addSrcRegisterForBranch(input_instr, srcCount, REG_STACK_POINTER, verbose, bType))
                return false;
            dstCount = 0;
            if (!addDstRegistersForBranch(input_instr, dstCount, {REG_STACK_POINTER, REG_INSTRUCTION_POINTER}, verbose, bType))
                return false;
            break;
        default:
            break;
    }
    return true;
}

// -----------------------------------------------------------------------------
// update_inst_registers: Update input_instr registers for the given instruction.
void update_inst_registers(void *dcontext, memref_t record, instr_t instr, input_instr &input_instr, bool verbose = false) {
    int srcCount = 0;
    int dstCount = 0;
    uint used_flag = instr_get_arith_flags(&instr, DR_QUERY_DEFAULT);

    for (int i = 0; i < instr_num_srcs(&instr); i++) {
        opnd_t opnd = instr_get_src(&instr, i);
        reg_id_t reg = opnd_get_reg_used(opnd, 0);
        if (verbose) {
            std::cout << "src register: " << reg << std::endl;
        }
        if (!addSrcRegister(input_instr, srcCount, reg, verbose))
            return;
    }

    for (int i = 0; i < instr_num_dsts(&instr); i++) {
        opnd_t opnd = instr_get_dst(&instr, i);
        reg_id_t reg = opnd_get_reg_used(opnd, 0);
        if (verbose) {
            std::cout << "dst register: " << reg << std::endl;
        }
        if (!addDstRegister(input_instr, dstCount, reg, verbose))
            return;
    }

    if (TESTANY(EFLAGS_WRITE_ARITH, used_flag)) {
        if (verbose) std::cout << "EFLAGS is written" << std::endl;
        if (!addDstRegisterForBranch(input_instr, dstCount, REG_FLAGS, verbose, BRANCH_INVALID))
            return;
    }
    if (TESTANY(EFLAGS_READ_ARITH, used_flag)) {
        if (verbose) std::cout << "EFLAGS is read" << std::endl;
        if (!addSrcRegisterForBranch(input_instr, srcCount, REG_FLAGS, verbose, BRANCH_INVALID))
            return;
    }

    if (record.instr.type == TRACE_TYPE_INSTR)
        return;

    {
        BranchType bType = getBranchType(record.instr.type);
        dstCount = 0;
        if (!addDstRegisterForBranch(input_instr, dstCount, REG_INSTRUCTION_POINTER, verbose, bType))
            return;
    }

    if (!assignBranchRegisters(input_instr, srcCount, dstCount, record.instr.type, verbose))
        return;
}

// -----------------------------------------------------------------------------
// update_branch_info: Update branch flags and fault counts.
void update_branch_info(memref_t record, input_instr &input_instr) {
    switch (record.instr.type) {
        case TRACE_TYPE_INSTR:
            input_instr.is_branch = false;
            input_instr.branch_taken = false;
            break;
        case TRACE_TYPE_INSTR_DIRECT_JUMP:
        case TRACE_TYPE_INSTR_INDIRECT_JUMP:
        case TRACE_TYPE_INSTR_DIRECT_CALL:
        case TRACE_TYPE_INSTR_INDIRECT_CALL:
        case TRACE_TYPE_INSTR_RETURN:
        case TRACE_TYPE_INSTR_TAKEN_JUMP:
            input_instr.is_branch = true;
            input_instr.branch_taken = true;
            break;
        case TRACE_TYPE_INSTR_UNTAKEN_JUMP:
            input_instr.is_branch = true;
            input_instr.branch_taken = false;
            break;
        case TRACE_TYPE_INSTR_CONDITIONAL_JUMP:
            input_instr.is_branch = true;
            input_instr.branch_taken = false;
            std::cerr << "unknown conditional jump" << std::endl;
            failure_counts[FAULT_UNKNOWN_CONDITIONAL_BRANCH_TAKEN]++;
            break;
        default:
            input_instr.is_branch = false;
            input_instr.branch_taken = false;
            std::cerr << "Unknown instruction type: " << record.instr.type << std::endl;
            failure_counts[FAULT_UNKNOWN_INST]++;
            break;
    }
}

// -----------------------------------------------------------------------------
// simulate_core: Processes records from a given stream.
void simulate_core(void* dcontext, scheduler_t::stream_t *stream, int thread_id, bool verbose,
                   SimStats &stats, const std::string &output_file_path,
                   const std::string &output_file_name) {

    std::unordered_map<size_t, instr_t> pc_instr_map;

    auto start_time = std::chrono::high_resolution_clock::now();
    std::string filename = output_file_path;
    if (!filename.empty() && filename.back() != '/' && filename.back() != '\\')
        filename += "/";

	int width = 4; // number of digits you want, e.g. 0001, 0002, ...
	std::ostringstream oss;
	oss << std::setw(width) << std::setfill('0') << thread_id;

    filename += output_file_name + "_" + oss.str() + ".champsim.gz";
    
    gzFile gz_out = gzopen(filename.c_str(), "w");
    if (!gz_out) {
        std::cerr << "Failed to open gz file " << filename << "\n";
        return;
    }
    
    if (verbose) {
        std::cout << "Thread " << thread_id << " processing filetype " 
                  << stream->get_filetype() << "\n";
    }
    
    auto decode_cache_ = init_decode_cache_and_set_flags(dcontext, stream);
    
    memref_t record;
    scheduler_t::stream_status_t status = stream->next_record(record);
    size_t local_count = 0;
    int64_t last_tid = -1245632134;
    
    while (status != scheduler_t::STATUS_EOF) {
        if (status == scheduler_t::STATUS_WAIT || status == scheduler_t::STATUS_IDLE) {
            std::this_thread::yield();
            status = stream->next_record(record);
            continue;
        }
        assert(status == scheduler_t::STATUS_OK);
        
        if (type_is_instr(record.instr.type)) {
            //size_t count = g_global_inst_count.fetch_add(1) + 1;
            //if (g_target_inst_count != 0 && count >= g_target_inst_count)
                //break;
	    local_count++;
            if (g_target_inst_count != 0 && local_count >= g_target_inst_count)
                break;
            //if (count % 100000 == 0) {
            if (local_count % 1000000 == 0) {
                auto elapsed_time = std::chrono::high_resolution_clock::now() - start_time;
                double elapsed_seconds = std::chrono::duration<double>(elapsed_time).count();
		print_mutex.lock();
                std::cout << "Processed " << local_count << " instructions, " << " Core " <<  thread_id << " cache size: " 
                          //<< decode_cache_->decode_cache_.size()
			  << "TODO"
			  << " Elapsed time: "
                          << elapsed_seconds << " seconds" << ", MI/s: "
                          << (local_count / 1000000.0) / elapsed_seconds << "\n";
		print_mutex.unlock();
            }
            stats.total_insts++;
            // if (decode_cache_->decode_cache_.size() > 200000) {
            //     std::cerr << "Cache size exceeded 200000, clearing cache\n";
            //     decode_cache_->decode_cache_.clear();
            // }

            if (verbose) {
                std::cout << std::left << std::setw(12) << "ifetch" << std::right
                          << std::setw(2) << record.instr.size << " byte(s) @ 0x" << std::hex
                          << std::setfill('0') << std::setw(sizeof(void *) * 2) << record.instr.addr
                          << std::dec << std::setfill(' ') << "\n";
                // std::cout << disasm_info->disasm_;
            }

	    int64_t cur_tid = record.instr.tid;
	    if (last_tid != cur_tid) {
		    // core id, instr count, thread id
		print_mutex.lock();
		std::cout << "***" << thread_id << "," << local_count << "," << cur_tid << "\n";
		print_mutex.unlock();
		last_tid = cur_tid;
	    }

            size_t pc = record.instr.addr;
            auto it = pc_instr_map.find(pc);
            instr_t instr;
            if (it != pc_instr_map.end()) {
                instr = it->second;
            } else {
                // disasm_info_t *disasm_info;
                // std::string error_string = decode_cache_->add_decode_info(record.instr, disasm_info);
                // instr = *disasm_info->instr;
                // pc_instr_map[pc] = instr;

                instr_init(dcontext, &instr);
                const app_pc decode_pc = reinterpret_cast<app_pc>(pc);
                app_pc nextpc = decode_from_copy(dcontext, record.instr.encoding, decode_pc, &instr);
                pc_instr_map[pc] = instr;
            }

            input_instr inst;
            update_inst_registers(dcontext, record, instr, inst, verbose);
            inst.ip = record.instr.addr;
            assert(inst.ip != 0);
            update_branch_info(record, inst);

            if (inst.is_branch) {
                stats.branch_count++;
                BranchType bType = getBranchType(record.instr.type);
                if (bType != BRANCH_INVALID)
                    stats.branch_type_counts[bType]++;
            }
            
            memref_t new_record;
            status = stream->next_record(new_record);
            int inst_source_memory = 0;
            int inst_dest_memory = 0;
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
                    inst.source_memory[inst_source_memory++] = new_record.instr.addr;
                    assert(new_record.instr.addr != 0);
                } else {
                    inst.destination_memory[inst_dest_memory++] = new_record.instr.addr;
                    assert(new_record.instr.addr != 0);
                }
                status = stream->next_record(new_record);
            }
            gzwrite(gz_out, &inst, sizeof(inst));
            record = new_record;
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
                   int num_cores) {
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
        !TESTANY(build_target_arch_type(), filetype_)) {
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

    all_stats.resize(num_cores);
    std::vector<std::thread> threads;
    threads.reserve(num_cores);
    for (int i = 0; i < num_cores; ++i) {
        threads.emplace_back([dcontext, i, &scheduler, verbose, &all_stats, &output_file_path, &output_file_name]() {
            simulate_core(dcontext, scheduler.get_stream(i), i, verbose, all_stats[i],
                          output_file_path, output_file_name);
        });
    }
    for (std::thread &thread : threads)
        thread.join();
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
    bool verbose = true;
    size_t target_count = 0;
    int num_cores = 8; // Default number of cores

    int opt;
    int option_index = 0;
    static struct option long_options[] = {
        {"trace_folder",     required_argument, 0, 't'},
        {"output_file_path", required_argument, 0, 'p'},
        {"output_file_name", required_argument, 0, 'f'},  // Changed from -n to -f.
        {"num_cores",        required_argument, 0, 'n'},  // New option for number of cores.
        {"quiet",            no_argument,       0, 'q'},
        {"target_count",     required_argument, 0, 'c'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "t:p:f:n:qc:", long_options, &option_index)) != -1) {
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
            case 'q':
                verbose = false;
                break;
            case 'c':
                target_count = std::stoull(optarg);
                break;
            default:
                std::cerr << "Usage: " << argv[0]
                          << " --trace_folder <folder> [--output_file_path <path>] [--output_file_name <name>] [--num_cores <num>] [--quiet] [--target_count <num>]\n";
                return 1;
        }
    }

    if (trace_directory.empty()) {
        std::cerr << "Error: --trace_folder is required.\n";
        std::cerr << "Usage: " << argv[0]
                  << " --trace_folder <folder> [--output_file_path <path>] [--output_file_name <name>] [--num_cores <num>] [--quiet] [--target_count <num>]\n";
        return 1;
    }

    g_target_inst_count = target_count;

    setrlim();

    std::vector<SimStats> thread_stats;
    auto start = std::chrono::high_resolution_clock::now();
    run_scheduler(trace_directory, verbose, thread_stats, output_file_path, output_file_name, num_cores);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    int total_insts = 0;
    int total_branches = 0;
    std::array<int, 8> total_branch_types = {0, 0, 0, 0, 0, 0, 0, 0};
    for (const auto &s : thread_stats) {
        total_insts += s.total_insts;
        total_branches += s.branch_count;
        for (int i = 0; i < 8; i++) {
            total_branch_types[i] += s.branch_type_counts[i];
        }
    }
    
    std::cout << "\n=== Simulation Statistics ===\n";
    std::cout << "Total instructions processed: " << total_insts << "\n";
    std::cout << "Total branch instructions:    " << total_branches << "\n";
    
    std::cout << "\nBranch Types for each core:\n";
    for (int i = 0; i < thread_stats.size(); i++) {
        std::cout << "Thread " << i << ":\n";
        for (int j = 0; j < 8; j++) {
            std::cout << branch_type_names[j] << ": " << thread_stats[i].branch_type_counts[j] << "\n";
        }
    }

    std::cout << "\nFaults:\n";
    for (int i = 0; i < FAULT_MAX; ++i) {
        std::cout << fault_names[i] << ": " << failure_counts[i] << "\n";
    }
    
    std::cout << "\nBranch Destination Register Overflows:\n";
    for (int i = 0; i < 8; i++) {
        std::cout << branch_type_names[i] << ": " << branch_dst_overflow_counts[i] << "\n";
    }
    
    std::cout << "\nBranch Source Register Overflows:\n";
    for (int i = 0; i < 8; i++) {
        std::cout << branch_type_names[i] << ": " << branch_src_overflow_counts[i] << "\n";
    }
    
    for (int i = 0; i < thread_stats.size(); i++) {
        std::cout << "Thread " << i << " processed " << thread_stats[i].total_insts << " instructions\n";
    }
    
    std::cout << "\nTime taken: " << elapsed.count() << " seconds\n";
    std::cout << "Seconds per million instructions: " << elapsed.count() / (total_insts * 1.0 / 1e6) << "\n";

    return 0;
}
