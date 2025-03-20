
디스크 정보 확인 
```
lsblk
```

smartctl 설치
```
sudo apt install smartmontools
sudo apt install smartmontools nvme-cli
```


smartctl 사용
```
sudo smartctl -A /dev/nvme0n1 > before_rocksdb.txt
sudo nvme smart-log /dev/nvme0n1 > before_rocksdb.txt

./db_bench --benchmarks=fillrandom --wal_dir= --db=
sudo smartctl -A /dev/WAl_디바이스 > after_rocksdb.txt
```

----

smartctl 예제

                Smart Log for NVMe device: /dev/nvme0n1
                Temperature:        34 Celsius
                Available Spare:    100%
                Available Spare Threshold:  10%
                Percentage Used:    5%
                Data Units Written: 123456789


nand block 확인 (논리)
```
sudo nvme id-ns /dev/nvme0n1
```


RocksDB 실행 중 102,848 LBA가 추가로 쓰임
- 15350796120 - 15350693272 = 102848 (LBA : 512 bytes)
- 102848 * 512 = 52,659,456 bytes
- 52,659,456 / 1,000,000,000 ≈ 0.0527 GB (약 52.7 MB)
- 52.7 MB / 3 MB ≈ 17.57 P/E 사이클 (TLC NAND 기준 블록 크기: 3MB)


