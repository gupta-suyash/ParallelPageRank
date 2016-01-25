#!/bin/bash
awk -F: '{printf("%s\n",$2)}' $1 > inputgraph.txt
g++ serial_matvec.cpp -o serial_matvec
./serial_matvec inputgraph.txt >pagerank.result
