#!/bin/bash

# 创建或重置CSV文件
output_file="db_bench_results.csv"
echo "Run,Fillrandom_microsec_per_op,Fillrandom_MB_per_sec,Fillrandom_value,Readrandom_microsec_per_op,Readrandom_MB_per_sec,Readrandom_value" > $output_file

# 重复运行50次
for i in {1..50}
do
    echo "Running iteration $i..."
    
    # 运行db_bench命令并将输出保存到临时文件
    ./db_bench --benchmarks="fillrandom,readrandom,stats" --num=40000000 --mod=7 --bwise=1 > temp_output.txt
    
    # 提取fillrandom和readrandom的结果
    fillrandom_result=$(grep "fillrandom   :" temp_output.txt)
    readrandom_result=$(grep "readrandom   :" temp_output.txt)
    
    # 使用正则表达式提取相关数值
    fillrandom_microsec_per_op=$(echo $fillrandom_result | awk '{print $3}')
    fillrandom_MB_per_sec=$(echo $fillrandom_result | awk '{print $5}')
    fillrandom_value=$(echo $fillrandom_result | awk '{print $7}')

    readrandom_microsec_per_op=$(echo $readrandom_result | awk '{print $3}')
    readrandom_MB_per_sec=$(echo $readrandom_result | awk '{print $5}')
    readrandom_value=$(echo $readrandom_result | awk '{print $10}')

    # 将结果添加到CSV文件
    echo "$i,$fillrandom_microsec_per_op,$fillrandom_MB_per_sec,$fillrandom_value,$readrandom_microsec_per_op,$readrandom_MB_per_sec,$readrandom_value" >> $output_file
done

# 删除临时文件
rm temp_output.txt

echo "All runs completed. Results saved to $output_file."
