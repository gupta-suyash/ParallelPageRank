#!/bin/bash
awk -F: '{printf("%s\n",$2)}' $1 > inputgraph.txt
g++ -std=c++0x -pthread -o pthread_matvec pthread_matvec.cpp
./pthread_matvec inputgraph.txt $2 $3 >pagerank.result
tail -n +3 pagerank.result > tmp.txt
sort -g -k2,2 tmp.txt > sorted_nodes.txt
tail -n $4 sorted_nodes.txt > highest_ids.txt
tac highest_ids.txt
