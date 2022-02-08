C1o=gcc
MAJOR = 1
MINOR = 9
BUILD = $(shell date +"%y%m%d.%H%M")
VERSIONA = "\"$(MAJOR).$(MINOR).$(BUILD)\""
VERSION = "$(MAJOR).$(MINOR).$(BUILD)"
VERSIONR = "$(MAJOR).$(MINOR)"
CPPFLAGS = -DVERSION=$(VERSIONA)

CFLAGS = -fopenmp -std=c99


all:  build debug prof
build: analysis lmp2sqt.simple analysis_binary
analysis_binary: analysis_binary.c snapshot.c common.c 
	$(CC) $(CFLAGS) analysis_binary.c snapshot.c common.c  -lm -O3 -o $@.$(VERSION).out -DNDEBUG
analysis: analysis.c snapshot.c common.c 
	$(CC) $(CFLAGS) analysis.c snapshot.c common.c  -lm -O3 -o $@.$(VERSION).out -DNDEBUG
lmp2sqt.simple: lmp2sqt.c snapshot.c common.c 
	$(CC) $(CFLAGS) lmp2sqt.c snapshot.c common.c -lm -O3  -o $@.$(VERSION).out -DNDEBUG
debug: debug_analysis  debug_lmp2sqt.simple debug_analysis_binary
prof: prof_lmp2sqt.simple 
debug_analysis_binary:analysis_binary.c snapshot.c common.c
	$(CC) $(CFLAGS) analysis_binary.c snapshot.c common.c -lm -ggdb -o $@.$(VERSION).gdb
debug_analysis:analysis.c snapshot.c common.c
	$(CC) $(CFLAGS) analysis.c snapshot.c common.c -lm -ggdb -o $@.$(VERSION).gdb
prof_lmp2sqt.simple:lmp2sqt.c snapshot.c common.c
	$(CC) $(CFLAGS) lmp2sqt.c snapshot.c common.c -lm  -pg -o $@.$(VERSION).gdb
debug_lmp2sqt.simple:lmp2sqt.c snapshot.c common.c
	$(CC) $(CFLAGS) lmp2sqt.c snapshot.c common.c -lm  -ggdb -o $@.$(VERSION).gdb

install: analysis.out lmp2sqt.simple.out analysis_binary.out
	cp analysis.out lmp2sqt.simple.out analysis_binary.out /home/kh/bin/
clean: clean_out clean_gdb clean_gprof

clean_out: *.out
	rm -r *.out 
clean_gdb: debug*gdb
	rm -r debug*.gdb
clean_gprof: prof*gdb
	rm -r prof*.gdb
