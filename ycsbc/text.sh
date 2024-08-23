#!/bin/bash

python3 /home/eros/workspace-lsm/Bourbon/test.py YSCB start

# 결과를 저장할 디렉토리 설정
output_directory="/mnt/1tb/ycsb/ycsb-test"

# 디렉토리가 존재하지 않으면 생성
if [ ! -d "$output_directory" ]; then
    mkdir -p "$output_directory"
fi

# 현재 날짜와 시간을 파일 이름에 사용
current_time=$(date "+%Y%m%d-%H%M%S")
output_file="$output_directory/unifrom-whisckey-20M-$current_time.txt"

# 명령 실행 및 결과 저장
sudo ./ycsbc input/20M > "$output_file"

echo "결과가 $output_file 에 저장되었습니다."

#work

python3 /home/eros/workspace-lsm/Bourbon/test.py YSCB end
