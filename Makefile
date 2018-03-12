all:  build debug
build: analysis.out lmp2sqt.simple.out analysis_binary.out
analysis_binary.out: analysis_binary.c snapshot.c common.c 
	gcc analysis_binary.c snapshot.c common.c  -lm -O3 -o $@ -DNDEBUG
analysis.out: analysis.c snapshot.c common.c 
	gcc analysis.c snapshot.c common.c  -lm -O3 -o $@ -DNDEBUG
lmp2sqt.simple.out: lmp2sqt.c snapshot.c common.c 
	gcc lmp2sqt.c snapshot.c common.c -lm -O3 -fopenmp -o $@ -DNDEBUG
debug: analysis.gdb lmp2sqt.simple.gdb analysis_binary.gdb
analysis_binary.gdb:analysis_binary.c snapshot.c common.c
	gcc analysis_binary.c snapshot.c common.c -lm -ggdb -o $@
analysis.gdb:analysis.c snapshot.c common.c
	gcc analysis.c snapshot.c common.c -lm -ggdb -o $@
lmp2sqt.simple.gdb:lmp2sqt.c snapshot.c common.c
	gcc lmp2sqt.c snapshot.c common.c -lm -fopenmp -ggdb -o $@

install: analysis.out lmp2sqt.simple.out analysis_binary.out
	cp analysis.out lmp2sqt.simple.out analysis_binary.out /home/kh/bin/

