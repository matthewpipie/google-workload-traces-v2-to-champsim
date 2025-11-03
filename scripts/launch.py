import os
import re
import subprocess
import time
import resource
import os

trace_base_path = "/mnt/storage/traces/gtrace_v2/external-traces-v2/"
output_dir = "/mnt/storage/traces/gtrace_v2_champsim/"
num_instructions_per_core = "100000000"

# The program will crash if the open limit cap isn't increased.
# Get the current soft and hard limits for open files
soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)
print(f"Current open file limits: soft={soft_limit}, hard={hard_limit}")
# Increase them
new_soft_limit = 50000
resource.setrlimit(resource.RLIMIT_NOFILE, (new_soft_limit, hard_limit))

# List all folders in the trace_base_path
trace_folders = [f for f in os.listdir(trace_base_path) if os.path.isdir(os.path.join(trace_base_path, f))]
print("Trace folders:")
print(trace_folders)

trace_folder_full_path = [os.path.join(trace_base_path, f) for f in trace_folders]
print("Trace folders full path:")
print(trace_folder_full_path)

core_map = {}

# Build core_map by reading the *_info.textproto file for each trace
for trace in trace_folders:
    description_file = os.path.join(trace_base_path, trace + "/aux/info.textproto")
    with open(description_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            pattern = r"peak_live_core_count: (\d+)"
            match = re.search(pattern, line)
            if match:
                core_map[trace] = match.group(1)
                break

print("Core map:")
print(core_map)

# Ensure the output directory exists

if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

# Create individual directories for each trace inside the output directory
for trace in trace_folders:
    trace_output_dir = os.path.join(output_dir, trace)
    os.makedirs(trace_output_dir, exist_ok=True)
    print(f"Created directory for {trace} in {output_dir}")

# Launch processes asynchronously and keep track of them along with their log file handles
processes = []  # will store tuples of (trace, process, log_file)

for trace in trace_folders:
    trace_output_dir = os.path.join(output_dir, trace)
    # Build the command
    if trace not in core_map:
        print(f"ERROR: {trace} not found in core_map")
        exit(1)
    command = [
        "../build/converter",
        "-t", os.path.join(trace_base_path, trace + "/trace"),
        "-c", num_instructions_per_core,
        "-n", core_map[trace],
        "-p", trace_output_dir,
        "-f", trace
    ]
    print("Command:")
    print(command)
    
    # Set up the log file path and open the file outside the context manager so it remains open
    log_file_path = os.path.join(trace_output_dir, f"{trace}.log")
    if os.path.exists(log_file_path):
        print(f"Log file {log_file_path} already exists. skipping.")
        continue
    
    log_file = open(log_file_path, "w")
    
    # Launch the process asynchronously, redirecting stdout and stderr to the log file
    process = subprocess.Popen(command, stdout=log_file, stderr=subprocess.STDOUT)
    processes.append((trace, process, log_file))

# Monitor processes and print a message as each one finishes
while processes:
    # Iterate over a copy of the list to safely remove items
    for trace, process, log_file in processes.copy():
        retcode = process.poll()
        if retcode is not None:  # Process has finished
            print(f"Process for trace '{trace}' finished with exit code {retcode}.")
            if retcode != 0:
                print(f"*** Check the log for more info! A crash has occured.")
            log_file.close()  # Close the log file for this process
            processes.remove((trace, process, log_file))
    time.sleep(0.5)  # Wait briefly before checking again

print("All processes have finished.")
