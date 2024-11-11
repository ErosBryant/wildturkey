#!/bin/bash

python3 ./test.py testing start 

# Define the desired --num values in an array
# nums=(64000000)
nums=(20000000)


# Define various configurations
# memtable_size=(4)
max_file_size=(2 8 16)
number_of_runs=2
# bwise=(1 0)
lacd=(1 2 3 4 5 6 7 8 9 10 15)
file_error=(2 4 8 16 32)

# fb_w wiki_w book_w
workload=(osm_w fb_w wiki_w book_w)
# lac=(5)
mod=(10)
# file_error=(22)

current_time=$(date "+%Y%m%d-%H%M%S")
# Define output directories
# output_dir="/mnt/lac-sec/ad-wt-bour/bourbon&wt-last/bourbon/"
output_dir="/mnt/wildturkey/real_workload2/"

test_dir="/home/eros/workspace-lsm/wildturkey/build/"

total_experiment="/mnt/1tb/lac_experiment/"



# Create output directories if they do not exist
if [ ! -d "$output_dir" ]; then
   mkdir -p "$output_dir"
   mkdir -p "${output_dir}summary_results/"
fi

# if [ ! -d "$total_experiment" ]; then
#    mkdir -p "$total_experiment"
# fi

# Execute the db_bench command for each configuration and save results
for num in "${nums[@]}"; do
   for wkload in "${workload[@]}"; do
         for md in "${mod[@]}"; do

            for i in $(seq 1 $number_of_runs); do

               output_file="${output_dir}mod=${md}-wkload=${wkload}-num=${num}_${i}.csv"
               
               echo "Running db_bench with --num=$num -wkload=${wkload}  --mod=${md} " > "$output_file"

               # Run the benchmark
               # uni40,uniread,stats
               # osm_w,real_r,stats
               # fillrandom,readrandom
               # --lac=$lacd 
               # --bwise=$bw
               # --max_file_size=$max
               # --lsize=${max/2} 
               # --file_error=$err
               # f=$((max / 2)) 
               # --lsize=$f
               ${test_dir}/db_bench --benchmarks="${wkload},real_r,stats" --mod=$md --num=$num >> "$output_file"
               echo "-------------------------------------" >> "$output_file"


         
               sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'
            done
         
         done
   done
done


python3 ./test.py testing end