
디스크 정보 확인 
```
lsblk
```

smartctl 설치
```
sudo apt install smartmontools
```


smartctl 사용
```
sudo smartctl -A /dev/WAl_디바이스 > before_rocksdb.txt
./db_bench --benchmarks=fillrandom --wal_dir= --db=
sudo smartctl -A /dev/WAl_디바이스 > after_rocksdb.txt
```

----

smartctl 예제
- db_bench 실행 전 
- 
            ID# ATTRIBUTE_NAME          FLAG     VALUE WORST THRESH TYPE      UPDATED  WHEN_FAILED RAW_VALUE
            *** 177 Wear_Leveling_Count     0x0013   099   099   000    Pre-fail  Always       -       5****
            *** 241 Total_LBAs_Written      0x0032   099   099   000    Old_age   Always       -       15350693272***




- db_bench 실행 후


            *** 177 Wear_Leveling_Count     0x0013   099   099   000    Pre-fail  Always       -       5 ***
            *** 241 Total_LBAs_Written      0x0032   099   099   000    Old_age   Always       -       15350796120 ***


RocksDB 실행 중 102,848 LBA가 추가로 쓰임
- 15350796120 - 15350693272 = 102848 (LBA : 512 bytes)
- 102848 * 512 = 52,659,456 bytes
- 52,659,456 / 1,000,000,000 ≈ 0.0527 GB (약 52.7 MB)
- 52.7 MB / 3 MB ≈ 17.57 P/E 사이클 (TLC NAND 기준 블록 크기: 3MB)


