
echo 3 > /proc/sys/vm/drop_caches

0-1-(2 & 18)-15-(5 & 19)-(12-14)
---
- Timer 0: level read
- Timer 1: file open (load ib+fb)
- Timer 2: model lookup (precdition)
- Timer 3: load datablock 
<!-- - Timer 3: key search time in file - first search -->
- Timer 15: FilteredLookup time
- Timer 5: load chunk & locate key
- TImer 12: Value reading time
- Timer 14: value read from memtable or immtable
----
- 还没设置
- Timer 17: get model time 
- Timer 18: search indexblock
- Timer 19: search datablock 


---
- Timer 6: Total key search time given files
- Timer 13 is the total time, which is the time we report.
- Timer 4 is the total time for all get requests.
- Timer 10 is the total time for all put requests.
- Timer 7 is the total compaction time.
- 
- Timer 11 is the total file model learning time.
- Timer 8 is the total level learning time
- timer 9: Total fresh write time (db load)
- timer 16: time to compact memtable

