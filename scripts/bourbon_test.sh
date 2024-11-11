#!/bin/bash

python3 ./test.py testing start 

# Define the desired --num values in an array
# nums=(64000000)
nums=(10000000)


# Define various configurations
# memtable_size=(4)
max_file_size=(2 8 16)
number_of_runs=1
# bwise=(1 0)
lacd=(1 2 3 4 5 6 7 8 9 10 15)
file_error=(2 4 8 16 32)

# fb_w wiki_w book_w
workload=(osm_w)
# lac=(5)
mod=(7)
# file_error=(22)

current_time=$(date "+%Y%m%d-%H%M%S")
# Define output directories
# output_dir="/mnt/lac-sec/ad-wt-bour/bourbon&wt-last/bourbon/"
output_dir="/mnt/motivation/$current_time/lac-sec/"

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
   # for bw in "${bwise[@]}"; do
   # for lac in "${lacd[@]}"; do
      # for max in "${max_file_size[@]}"; do
      for err in "${file_error[@]}"; do
         for md in "${mod[@]}"; do
         # Initialize summary output file
            summary_output="${output_dir}summary_results/err=${err}lac=${lac}-mod=${md}-num=${num}.csv"
            echo "  num,  run, write_micros/op, read_micros/op, write_MB/s, read_MB/s, lac , err,  mod , waf, memtable_stall, L0_stall, L0_slow_stall, avg_segment_size" > "$summary_output"


            # 변수 리스트 초기화
            write_micros_list=()
            read_micros_list=()
            write_mb_list=()
            read_mb_list=()
            for i in $(seq 1 $number_of_runs); do
               # Define output file
               # lac=${lacd}
               # max_file_size=${max}
               # error=${err}
               # bwise=${bw}
               output_file="${output_dir}mod=${md}ERR=${err}_lac=${lac}-num=${num}_${i}.csv"
               
               echo "Running db_bench with --num=$num  lac=$lac max=$max --mod=${md} " > "$output_file"

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
               ${test_dir}/db_bench --benchmarks="osm_w,real_r,stats" --mod=$md --file_error=$err --num=$num >> "$output_file"
               echo "-------------------------------------" >> "$output_file"


               # Extract performance data
               write_micros_per_op=$(grep "${wkload}" "$output_file" | awk '{for(i=1;i<=NF;i++) if($i=="micros/op;") print $(i-1)}')
               read_micros_per_op=$(grep "real_r" "$output_file" | awk '{for(i=1;i<=NF;i++) if($i=="micros/op;") print $(i-1)}')
               write_mb_per_s=$(grep "${wkload}" "$output_file" | awk '{for(i=1;i<=NF;i++) if($i=="MB/s;") print $(i-1)}')
               read_mb_per_s=$(grep "real_r" "$output_file" | awk '{for(i=1;i<=NF;i++) if($i=="MB/s") print $(i-1)}')
               waf=$(grep 'waf:' "$output_file" | awk -F':' '{print $2}')
               memtable_stall=$(grep 'memtable stall time' "$output_file" | awk '{print $(NF-1)}')
               l0_stall=$(grep 'L0 stall time' "$output_file" | awk '{print $(NF-1)}')
               l0_slow_stall=$(grep 'L0 slow stall time' "$output_file" | awk '{print $(NF-1)}')
               avg_segment_size=$(grep 'Average Segement Size' "$output_file" | awk '{print $NF}')
               
                        # 리스트에 성능 데이터 추가
               write_micros_list+=($write_micros_per_op)
               read_micros_list+=($read_micros_per_op)
               write_mb_list+=($write_mb_per_s)
               read_mb_list+=($read_mb_per_s)
               # Append data to summary output file
               echo "$num, $i,     $write_micros_per_op,          $read_micros_per_op,         $write_mb_per_s,       $read_mb_per_s,  $lac , $err,  $mod   ,  $waf,     $memtable_stall,    $l0_stall,      $l0_slow_stall,          $avg_segment_size" >> "$summary_output"

               # Clear system cache
               sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'
         done

       
         avg_write_micros=$(echo "${write_micros_list[@]}" | awk '{sum=0; for(i=1;i<=NF;i++) sum+=$i; print sum/NF}')
         avg_read_micros=$(echo "${read_micros_list[@]}" | awk '{sum=0; for(i=1;i<=NF;i++) sum+=$i; print sum/NF}')
         avg_write_mb=$(echo "${write_mb_list[@]}" | awk '{sum=0; for(i=1;i<=NF;i++) sum+=$i; print sum/NF}')
         avg_read_mb=$(echo "${read_mb_list[@]}" | awk '{sum=0; for(i=1;i<=NF;i++) sum+=$i; print sum/NF}')

         # 평균값을 summary_output 파일에 추가
         echo "Average, avg_write_micros, avg_read_micros, avg_write_mb, avg_read_mb" >> "$summary_output"
         echo "Average, $avg_write_micros, $avg_read_micros, $avg_write_mb, $avg_read_mb" >> "$summary_output"
         done
      done
   done
done


python3 ./test.py testing end
