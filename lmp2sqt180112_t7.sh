#!/bin/bash -
#===============================================================================
#
#          FILE: lmp2sqt180112.sh
#
#         USAGE: ./lmp2sqt180112.sh
#
#   DESCRIPTION: 
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: KIM Hyeok (kh), ekh0324@gmail.com
#  ORGANIZATION: Konkuk University
#       CREATED: 2018년 01월 12일 18시 29분 30초
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error
sleep 1
export OMP_NUM_THREADS=7
t1=$(date -u +%s)
exec >all.out
/home/kh/bin/lmp2sqt180112.out $@
t2=$(date -u +%s)
TIMER=$(date -u  -d "@$(( t2 - t1 ))" +"%H:%M:%S"  )
echo "The total execution time of the programe is $TIMER" >>TIME
