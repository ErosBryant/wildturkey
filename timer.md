Timer 0: Finde file (level에서 해당 file 찾기)
Timer 1: load,  reading time file(?) table cache에서 file 찾기
Timer 2: file model 찾기
Timer 3: 0
Timer 4: 0
Timer 5: 0
Timer 6: 0
Timer 7: 0
Timer 8: 0
Timer 9: 0
Timer 10: 0
Timer 11: 0
Timer 12: 0
Timer 13: 0
Timer 14: 0
Timer 15: 0
Timer 16: 0
Timer 17: 0
Timer 18: 0
Timer 19: 0

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
