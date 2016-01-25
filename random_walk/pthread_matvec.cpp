#include <cmath>
#include <cstdio>
#include <vector>
#include<string>
#include<stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <iterator>
#include <pthread.h>
#include<sys/time.h>
#include <unistd.h>
#include <random>

using namespace std;

// Defining global data.
int partitions = 16;			// Determines the number of threads to run. 
int steps = 119;
long size=0, vr=0;//, part[partitions];
std::vector<int> col_inds, row_ptrs;
std::vector<long> result, vals;
pthread_barrier_t barrier;
pthread_mutex_t minimum_value_lock;

void read_graphs(char *argv[]) {
	long i,j;
	int x; 
	string line;

	// Reading the graph file.
	ifstream fp;
	//printf("%s ",argv[1]);
	fp.open(argv[1]);			// Name of the input file
	//cout<<argv[1];

	std::getline(fp,line);
	std::stringstream stream(line);

	while(1) {
		int n;
		stream >> n;
		if(!stream)
			break;
		vals.push_back(n);	
	}
	size = vals.size();

	std::getline(fp,line);
	std::stringstream stream1(line);

	while(1) {
		int n;
		stream1 >> n;
		if(!stream1)
			break;
		col_inds.push_back(n);	
	}

	std::getline(fp,line);
	std::stringstream stream2(line);

	while(1) {
		int n;
		stream2 >> n;
		if(!stream2)
			break;
		row_ptrs.push_back(n);	
	}
	vr = row_ptrs.size();

	fp.close();


//	printf("\n Input vals: ");
//	for(i=0; i<*size; i++) {
//		printf("%d  ",vals[i]);
//	}
//
//	printf("\n Input cols: ");
//	for(i=0; i<*size; i++) {
//		printf("%d  ",col_inds[i]);
//	}
//
//	printf("\n Input rows: ");
//	for(i=0; i<*vr; i++) {
//		printf("%d  ",row_ptrs[i]);
//	}
//	printf("\n");

	//Generating Vector.	//printf("Enter vector: \n");
	long sum=0;
	for(i=0; i<vr-1; i++) {
		result.push_back(0);
		sum++;
	}

	partitions 	= atoi(argv[2]);	// Getting the number of threads.
	steps		= atoi(argv[3]);	// Getting the steps.
}




void *perform_randomwalk(void *arg) {
	long myst = ((long*)arg)[1], myend = ((long*)arg)[2], myid = ((long*)arg)[0];
	long mysz = myend-myst;
	std::vector<long> lwalk;
	int i=0, j=0; 
	int currentNode, neighbors, nextNode, targetColumn,rowid;

	lwalk.resize(vr-1);

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0,vr-1);
	for(i=0; i<vr-1; i++)
		lwalk[i] = 0;

//	cout << "\n " << myid << " : start: " << myst << " : end: " << myend; 
//	fflush(stdout);
	
	struct timeval start, end;
	pthread_barrier_wait (&barrier);

	if(myid == 0) {
		gettimeofday(&start, NULL);		
	}

	for(i=myst; i<myend; i++) {
		currentNode = i;			// Starting node.
		for(j=0; j<steps; j++) {
			rowid		= row_ptrs[currentNode];
			neighbors 	= row_ptrs[currentNode+1]-rowid; 
			nextNode 	= int(distribution(generator))%neighbors;//rand()%neighbors;
			targetColumn 	= col_inds[rowid+nextNode];
			lwalk[targetColumn] += 1;
			currentNode	= targetColumn;	
		}
	}

	//pthread_barrier_wait (&barrier);

	while(1) {
		if(pthread_mutex_lock(&minimum_value_lock) == 0) {
			//cout << "\n I got lock: " << myid;
			for(i=0; i<vr-1; i++) {
				result[i] += lwalk[i];
			}

			//cout << "\n I am leaving the lock: " << myid;
			//fflush(stdout);

			pthread_mutex_unlock(&minimum_value_lock);
			break;
		}
	}
	
	pthread_barrier_wait (&barrier);

	if(myid == 0) {
		//cout << k << "\n";

		gettimeofday(&end, NULL);
		printf("time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec)
		  - (start.tv_sec * 1000000 + start.tv_usec)));
		for(i=0; i<vr-1; i++)
			cout << "\n" << i << " " << result[i];
	}
}


int main(int argc, char *argv[]) {
	int i;	long j;

	// Reading graphs.
	read_graphs(argv);

//	cout << "\n Vals: "; //<< vals.size();
//	for(i=0; i<size; i++) {
//		printf("%f ",nvals[i]);
//	}
//	
//	cout << "\n Cols: ";
//	for(i=0; i<vr; i++) {
//		printf("%d ",ncols[i]);
//	}
//	
//	cout << "\n Rows: ";
//	for(i=0; i<size; i++) {
//		printf("%d ",nrows[i]);
//	}

	// Creating partitions.	
	int part[partitions];

	//cout << "Steps: " << steps << " Threads: " << partitions;

	for(i=0; i<partitions; i++) 
		part[i] = 0;
	j = (vr-1)%partitions;
	for(i=1; i<partitions; i++) {
		part[i] = part[i-1] + (vr-1)/partitions;
	}

	pthread_t p_threads[partitions]; 
	pthread_barrier_init(&barrier, NULL, partitions);	
	pthread_mutex_init(&minimum_value_lock, NULL);

	for (i=0; i<partitions; i++) {
		long *arg_send = new long(3);
		arg_send[0] = i;			// Thread id.
		arg_send[1] = part[i];			// Thread index start.
		if(i != partitions-1 ) {		// Thread index end.
			arg_send[2] = part[i+1];
		}
		else {
			arg_send[2] = part[i]+((vr-1)/partitions)+j;	// All the remaining elements.
		}

		pthread_create(&p_threads[i], NULL, perform_randomwalk, (void *)arg_send);
	}

	for (i=0; i<partitions; i++) {
		pthread_join(p_threads[i], NULL);
	}
	
//	// Printing Result.
//	printf("\n Result: \n");
//	for(i=0; i<vr-1; i++) {
//		printf("%d ",result[i]);
//	}
//
//	printf("\n Bye \n");
}
