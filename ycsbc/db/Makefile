CC=g++
CFLAGS=-std=c++11 -c -g -Wall -fgnu-tm


INCLUDES=-I../ -I/home/eros/workspace-lsm/wildturkey/include/ -I/home/eros/workspace-lsm/wildturkey/mod/ -I/home/eros/workspace-lsm/wildturkey/mod/




SOURCES=$(wildcard *.cc core/*.cc db/*.cc)

OBJECTS=$(SOURCES:.cc=.o)


all: $(SOURCES) $(OBJECTS)

.cc.o:
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) $(OBJECTS)