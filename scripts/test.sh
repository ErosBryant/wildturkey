

# Define the desired --num values in an array
nums=(10000000)

# error_bound=(2)
# 2 4 8 16 32

number_of_runs=3

value_size=(100 500)
# sst_size=(1 2 3 4 8 16 32)
# sst_size=(2)

# 1 2 3 
sst_base_size=(1 2 3 4 8 16)

mods=(5)

workload=(fb_w)

# use_bwises=(0)

# Define output directories
output_dir="/home/eros/workspace-lsm/wildturkey/0909-leveldb/"

test_dir="/home/eros/workspace-lsm/wildturkey/build/"

current_time=$(date "+%Y%m%d-%H%M%S")

db_data_dir="/home/eros/workspace-lsm/wildturkey/db_dir/"

# Create output directories if they do not exist
if [ ! -d "$output_dir" ]; then
   mkdir -p "$output_dir"
fi

if [ ! -d "$db_data_dir" ]; then
   mkdir -p "$db_data_dir"
fi


# Execute the db_bench command for each configuration and save results
for num in "${nums[@]}"; do
   for mod in "${mods[@]}"; do
      # for error in "${error_bound[@]}"; do
      # for sst in "${sst_size[@]}"; do
      # for rwork in "${workload[@]}"; do
      for value in "${value_size[@]}"; do
         for sst_base in "${sst_base_size[@]}"; do
         # for bwises in "${use_bwises[@]}"; do
            for i in $(seq 1 $number_of_runs); do
               # output_file="${output_dir}SST_${num}_${mod}_--bwise=${bwises}_${rwork}_${i}_.txt"
               output_file="${output_dir}_${sst_base}SST_${num}_${mod}_${i}.txt"
               # echo "Running db_bench with --num=$num --sst_size=$sst --error_bounds=$error --times=$i" >> "$output_file"
               echo "Running db_bench with --num=$num  --mod=$mod --times=$i --sst_base=${sst_base} --value_size=${value}" >> "$output_file"
               ${test_dir}/db_bench --benchmarks="fillrandom,readrandom,stats" --max_file_size=${sst_base}  --value_size=$value --mod=$mod --num=$num --db=$db_data_dir >> "$output_file"
               # ${test_dir}/db_bench --benchmarks="${rwork},real_r,stats"  --max_file_size=${sst_base} --mod=$mod --num=$num --db=$db_data_dir >> "$output_file"
               echo "-------------------------------------" >> "$output_file"
            
               sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'
            done
         done
         
      done         
   done
done
