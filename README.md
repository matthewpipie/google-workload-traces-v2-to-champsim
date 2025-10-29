# google-workload-traces-v2-to-champsim
A conversion tool to generate Champsim traces from Google's Workload Traces (Version 2)

## Building

First, build dynamorio:
`cd dynamorio && mkdir build && cd build && cmake .. && make -j && cd ../..`

Then, build the conversion script:
`mkdir build && cd build && cmake .. && make -j`

You should see the output program `converter`. 

There is no need to build ChampSim, as we just need its header files.

## Running

# TODO
* type safe rand
* fix up rand
* fix up entire schedule function
* create per-thread args
* turn quiet into verbose
* download google trace instructions
* usage instructions