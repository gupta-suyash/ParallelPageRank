#!/bin/bash
awk -F: '{printf("%s\n",$2)}' $1 > inputgraph.txt
gcc serial_matvec.c -lm -o serial_matvec
./serial_matvec inputgraph.txt >pagerank.result
