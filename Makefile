build:
	gcc analysis.c snapshot.c common.c  -lm -O3 -o analysis.out
debug:
	gcc analysis.c snapshot.c common.c -lm -ggdb -o analysis.gdb
