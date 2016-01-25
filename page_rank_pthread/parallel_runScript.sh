#!/bin/bash
awk -F: '{printf("%s\n",$2)}' $1 > inputgraph.txt
g++ pthread_matvec_working.cpp -pthread -o pthread_matvec_working
./pthread_matvec_working inputgraph.txt >pagerank.result
