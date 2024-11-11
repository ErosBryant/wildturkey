CC=/usr/bin/g++
# CFLAGS=-std=c++11 -g -Wall -pthread -I./
CFLAGS=-std=c++11 -g -Wall -pthread -I./ -I/home/eros/workspace-lsm/wildturkey/include/  -I/home/eros/workspace-lsm/wildturkey/mod/

# -lleveldb
# LDFLAGS= -lpthread  -lz -lsnappy
# LDFLAGS= -lpthread -L/home/eros/workspace-lsm/Bourbon/build -lz -lsnappy
LDFLAGS= -lpthread /home/eros/workspace-lsm/wildturkey/build/libleveldb.a -lz -lsnappy


SUBDIRS= core db 
SUBSRCS=$(wildcard core/*.cc) $(wildcard db/*.cc)
OBJECTS=$(SUBSRCS:.cc=.o)
EXEC=ycsbc

all: $(SUBDIRS) $(EXEC)

$(SUBDIRS):
	$(MAKE) -C $@

$(EXEC): $(wildcard *.cc) $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir $@; \
	done
	$(RM) $(EXEC)

.PHONY: $(SUBDIRS) $(EXEC)
