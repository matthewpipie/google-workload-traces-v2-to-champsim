# google-workload-traces-v2-to-champsim
A conversion tool to generate Champsim traces from Google's Workload Traces (Version 2)

## Downloading Google Traces

Follow instructions on https://dynamorio.org/google_workload_traces.html#sec_google_get . 

You may find this command useful for recursively downloading all traces: https://docs.cloud.google.com/storage/docs/downloading-objects#cli-download-object .

The traces are approximately 842GB.

## Building the Converter

First, build DynamoRIO:
`cd dynamorio && mkdir build && cd build && cmake .. && make -j && cd ../..`

If you get errors here, be sure to init all submodules in this repository and in the DynamoRIO repository.

Then, build the conversion tool:
`mkdir build && cd build && cmake .. && make -j`

You should see the output program `converter`. 

There is no need to build ChampSim, as we just need its header files.

## Running

`./converter` accepts a few args:
* `--trace_folder` defines the path to the Google workload folder, e.g. `/mnt/storage/traces/gtrace_v2/external-traces-v2/arizona/trace`
* `--output_file_path` defines the folder where output files will be placed
* `--output_file_name` defines a prefix for the output files. I normally use the workload name
* `--num_cores` defines the number of traces to output. Each output trace is one "core" that has been reconstructed from Google's thread level traces
* `--verbose` enables debug-level output
* `--percore_target_count` defines the per-core cap on instructions converted before quitting. Set to 0 to ignore.
* `--global_target_count` defines the global cap on instructions converted before quitting. Set to 0 to ignore.

Notably, Google recommends simulating their threads using `peak_live_core_count` cores. So, you should set `--num_cores` to the `peak_live_core_count` inside `external-traces-v2/<workload>/aux/info.textproto`.

As an example, one might write:
`./converter --trace_folder /mnt/storage/traces/gtrace_v2/external-traces-v2/arizona/trace --output_file_path /mnt/storage/traces/gtrace_v2_champsim/arizona --output_file_name arizona --num_cores 6 --percore_target_count 1000000000 > /mnt/storage/traces/gtrace_v2_champsim/arizona/conversion.log`
in order to convert 1B instructions per-core on the arizona workload.

If you want to convert threads to their own ChampSim thread-level trace, you could set the `--trace_folder` to the actual thread file.

## Scripting

We provide a script, `launch.py`, which converts all of the workloads in parallel according to `peak_live_core_count`. You can find it in `scripts/`. Edit the paths accordingly before running, and the target instruction count, defined at the top of the file.

# Disclaimers and other notes

Clang is not supported by DynamoRIO, you should run these tools on GNU Linux with g++, I've never tried them on anything else and they probably don't work.

There are several sources of inaccuracies between the real Google workloads and the Google traces.
* The addition of tracing artificially increases the frequency of context switches above what may normally happen.
* The obsfucation preformed by Google combines several registers and anonymizes them, including the stack register, making representative simulations hard. More on this in a big comment near the top of main.cpp.
* More that I can't think of right now.

There are several sources of inaccuracies between the Google traces and the outputted ChampSim traces.
* Sometimes, instructions have too many source or destination registers or memory references. This occurs due to Google's obsfucation. To make the simulations slightly more accurate, we could increase ChampSim's `NUM_INSTR_SOURCES` and `NUM_INSTR_DESTINATIONS` in `trace_instruction.h`. However, we'd need our own fork of ChampSim to work with these traces, then.
* ChampSim fully depends on register dependencies in order to construct accurate timings between memory operations. Since Google has obfuscated these with combinations, our timings aren't going to be very accurate.
* The DynamoRIO scheduler expects traces to be parsed in conjunction with ChampSim, i.e. ChampSim should integrate with DynamoRIO's libraries to get the next instruction. This allows DynamoRIO access to the clock, allowing it to more accurately reconstruct the thread ordering across cores. However, this isn't possible in practice, because ChampSim cannot support multithreaded workloads. In this converter, we use DynamoRIO's default mode of operation, which assumes an 2.0 GHz processor with 0.5 IPC and a time quantum of 5 ms between preemptions (unless the thread gets blocked) (source: scheduler.h:800). This makes the outputted core traces not perfectly line up with their real-world counterpart, because time runs at a slightly different rate in this scheduler.
* More that I can't think of right now.

There are several sources of inaccuracies between ChampSim traces and ChampSim output:
* ChampSim only supports single-threaded simulations. Although newer ChampSim supports multicore, these cores are assumed to not interfere / to belong to different processes. The Google cores all run code in the same multithreaded process. Thus, the single-threaded simulations do not accurately represent sharing behavior between cores, or contention that would be present in the shared microarchitecture (e.g. LLC, DRAM).
* ChampSim isn't accurate in many ways.
* More that I can't think of right now.

All that said, I think this is roughly the best we can do right now. Perhaps integrating with gem5 or other simulators could be more fruitful, but there are already plenty of non-ChampSim causes of inaccuracies that I don't know how different it would be. Definitely worth a shot at some point!

There are some weird outputs in ChampSim with these traces, namely:
"[BTB] WARNING: target of return is a lower address than the corresponding call. This is usually a problem with your trace."
I'm still looking into this.

Each core tends to behave fairly differently, i.e. threads are not randomly multiplexed onto cores. Sampling only a subset of the cores to simulate in ChampSim is a bad idea, as is scheduling onto a different number of cores than Google recommends.

The expected output IPC is low: less than 1 on every workload except `arizona`. Good luck. 

# Credits

The original converter was written by Kan Zhu. It has been heavily reworked by Matthew Giordano, who now maintains this repo.

Thanks to Akanksha Jain and Enrico Deiana from Google for their help in writing this code.