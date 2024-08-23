#!/bin/bash

python3 ../test.py testing start 

# Define the desired --num values in an array
nums=(80000000)

# Define error bounds
error_bound=(8)

# Define various configurations
memtable_size=(4 8 16 32 64)
max_file_size=(2 4 8 16 32)
sst_time=(1)
number_of_runs=1
mod=(7)

# Define output directories
output_dir="/mnt/1tb/sec5_5_bour/"

test_dir="/home/eros/workspace-lsm/wildturkey/build/db_bench/"

total_experiment="/mnt/1tb/result/sec5_5_bour/"

current_time=$(date "+%Y%m%d-%H%M%S")

# Create output directories if they do not exist
if [ ! -d "$output_dir" ]; then
   mkdir -p "$output_dir"
fi

if [ ! -d "$total_experiment" ]; then
   mkdir -p "$total_experiment"
fi

# Execute the db_bench command for each configuration and save results
for num in "${nums[@]}"; do
   for sst_times in "${sst_time[@]}"; do
      for error in "${error_bound[@]}"; do
         for mem_size in "${memtable_size[@]}"; do
            for sst_size in "${max_file_size[@]}"; do
               declare -a fill_array=()
               declare -a fread_array=()
               declare -a segment_array=()

               # Initialize average variables
               fill_avg=0
               read_avg=0
               segment_avg=0

               # Initialize summary output file
               summary_output="${output_dir}summary_results_SST_${sst_times}.txt"
               echo "num, error_bound, run, fillrandom, readrandom, sst_times, segment_size" >> "$summary_output"

               for i in $(seq 1 $number_of_runs); do
                  output_file="${output_dir}SST_MEM_${mem_size}MB_${sst_size}MB_${sst_times}_Entry_${num}_error_${error}_${i}.txt"
                  echo "Running db_bench with --num=$num --error_bound=$error --sst_times=$sst_times" >> "$output_file"
                  ${test_dir}/db_bench --benchmarks="fillrandom,readrandom,stats" --mod=7 --num=$num --write_buffer_size=$mem_size --max_file_size=$sst_size --file_error=$error --sst_times=$sst_times >> "$output_file"
                  echo "-------------------------------------" >> "$output_file"
               
                  # Extract fillrandom and readrandom performance data
                  fillrandom=$(grep "fillrandom" "$output_file" | tail -1 | awk '{print $3}')
                  readrandom=$(grep "readrandom" "$output_file" | tail -1 | awk '{print $3}')
                  segment_size=$(grep "segement_size" "$output_file" | awk '{print $2}')

                  fill_array+=($fillrandom)
                  fread_array+=($readrandom)
                  segment_array+=($segment_size)
               
                  echo "$num, $error, $i, $fillrandom, $readrandom, $sst_times, $segment_size" >> "$summary_output"

                  sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'
               done
         
               # Calculate averages and write to file
               fill_avg=$(IFS=+; bc <<< "scale=2;(${fill_array[*]})/$number_of_runs")
               read_avg=$(IFS=+; bc <<< "scale=2;(${fread_array[*]})/$number_of_runs")
               segment_avg=$(IFS=+; bc <<< "scale=2;(${segment_array[*]})/$number_of_runs")
            
               echo "fillrandom_avg: $fill_avg" >> "$summary_output"
               echo "readrandom_avg: $read_avg" >> "$summary_output"
               echo "segment_avg: $segment_avg" >> "$summary_output"

               fill_avg_all=$(grep "fillrandom_avg" "$summary_output" | awk '{print $2}')
               read_avg_all=$(grep "readrandom_avg" "$summary_output" | awk '{print $2}')
               segment_avg_all=$(grep "segment_avg" "$summary_output" | awk '{print $2}')

               echo "SST: ${sst_times} _ error:${error}" >>  "${total_experiment}${current_time}_result.txt"
               echo "fillrandom_avg_all: $fill_avg_all" >> "${total_experiment}${current_time}_result.txt"
               echo "readrandom_avg_all: $read_avg_all" >> "${total_experiment}${current_time}_result.txt"
               echo "segment_avg_all: $segment_avg_all" >>"${total_experiment}${current_time}_result.txt"
      
            done
         done
      done
   done
done

python3 ../test.py testing end
