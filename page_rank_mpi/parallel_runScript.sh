#!/bin/bash
awk -F: '{printf("%s\n",$2)}' $1 > inputgraph.txt
mpicc MPI_matvec_new.c -o MPI_matvec_new
mpirun -machinefile machines -np $3 MPI_matvec_new inputgraph.txt $2 $3 >pagerank.result
