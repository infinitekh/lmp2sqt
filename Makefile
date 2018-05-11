CC=gcc
MAJOR = 1
MINOR = 0
BUILD = $(shell date +"%y%m%d.%H%M")
VERSIONA = "\"$(MAJOR).$(MINOR).$(BUILD)\""
VERSION = "$(MAJOR).$(MINOR).$(BUILD)"
CPPFLAGS = -DVERSION=$(VERSIONA)




all:  build debug
build: analysis.out lmp2sqt.simple.out analysis_binary.out
analysis_binary.out: analysis_binary.c snapshot.c common.c 
	$(CC) $(CFLAGS) analysis_binary.c snapshot.c common.c  -lm -O3 -o $@.$(VERSION) -DNDEBUG
analysis.out: analysis.c snapshot.c common.c 
	$(CC) $(CFLAGS) analysis.c snapshot.c common.c  -lm -O3 -o $@.$(VERSION) -DNDEBUG
lmp2sqt.simple.out: lmp2sqt.c snapshot.c common.c 
	$(CC) $(CFLAGS) lmp2sqt.c snapshot.c common.c -lm -O3 -fopenmp -o $@.$(VERSION) -DNDEBUG
debug: analysis.gdb lmp2sqt.simple.gdb analysis_binary.gdb
analysis_binary.gdb:analysis_binary.c snapshot.c common.c
	$(CC) $(CFLAGS) analysis_binary.c snapshot.c common.c -lm -ggdb -o $@.$(VERSION)
analysis.gdb:analysis.c snapshot.c common.c
	$(CC) $(CFLAGS) analysis.c snapshot.c common.c -lm -ggdb -o $@.$(VERSION)
lmp2sqt.simple.gdb:lmp2sqt.c snapshot.c common.c
	$(CC) $(CFLAGS) lmp2sqt.c snapshot.c common.c -lm -fopenmp -ggdb -o $@.$(VERSION)

install: analysis.out lmp2sqt.simple.out analysis_binary.out
	cp analysis.out lmp2sqt.simple.out analysis_binary.out /home/kh/bin/

