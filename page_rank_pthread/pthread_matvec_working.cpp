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

using namespace std;

// Defining global data.
const int partitions = 128;			// Determines the number of threads to run. 
int size, vr, part[partitions];
double total_sum[partitions];
std::vector<int> col_inds, row_ptrs;
std::vector<double> vec, vals;
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
		vec.push_back(1);
		sum++;
	}

	for(i=0; i<vr-1; i++) 
		vec[i] = vec[i]/sum;

//	cout << "\n Vector: ";
//	for(i=0; i<*vr-1; i++) {
//		cout << " " << vec[i];
//	}
}


void *perform_matvec(void *arg) {
	long i,j,flag=1, k=0, lindex;
	long myst = ((int*)arg)[1], myend = ((int*)arg)[2], myid = ((int*)arg)[0];
	long mysz = myend-myst;
	double temp, sum, check = 0.00001;

	std::vector<double> result;
        result.resize(mysz);

	for(i=0; i<mysz; i++) 
		result[i] = 0;

	struct timeval start, end;
	if(myid == 0) {
		gettimeofday(&start, NULL);
	}

	while(flag == 1) {
		sum = 0;	temp=0;
		lindex=0;
		for(i=myst; i<myend; i++) {
			temp = 0;
			for(j=row_ptrs[i]; j<row_ptrs[i+1]; j++) 
				temp += vals[j] * vec[col_inds[j]];
			result[lindex] = temp;
			lindex++;
		}

		lindex=0;
		for(i=myst; i<myend; i++) {
			temp = result[lindex] - vec[i];
			temp = pow(temp,2);
			sum += temp;
			lindex++;
		}

		total_sum[myid] = sum;

		// Barrier to let everyone to write sum.
		pthread_barrier_wait (&barrier);

		temp = 0;
		for(i=0; i<partitions; i++) {
			temp += total_sum[i];
		}

		sum = sqrt(temp);
		if(sum < check) {
			flag = 0;
		}

		k++;

	//	lindex = 0;
	//	for(i=myst; i<myend; i++) {
	//		vec[i] = result[lindex];
	//		lindex++;
	//	}

	//	if(myid == 0 || myid == 2 || myid == 4 || myid == 6 || myid == 8 || myid == 10 || myid == 12 ||
	//	myid == 14 || myid == 16 || myid == 18 || myid == 20 || myid == 22 || myid == 24 || myid == 26 ||
	//	myid == 28 || myid == 30 //|| myid == 18 || myid == 20 || myid == 22 || myid == 24 || myid == 26 ||
		if(myid % 2 == 0
		) {
			lindex=0;
			for(i=myst; i<myend; i++) {
				vec[i] = result[lindex];
				lindex++;
			}
		}

	//	else if(myid == 2) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 4) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 6) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 8) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 10) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 12) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 14) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 16) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 18) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 20) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 22) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 24) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 26) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 28) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 30) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

		pthread_barrier_wait (&barrier);

	//	if(myid == 1 || myid == 3 || myid == 5 || myid == 7 || myid == 9 || myid == 11 || myid == 13 ||
	//	myid == 15 || myid == 17 || myid == 19 || myid == 21 || myid == 23 || myid == 25 || myid == 27 ||
	//	myid == 29 || myid == 31 //|| myid == 18 || myid == 20 || myid == 22 || myid == 24 || myid == 26 ||
		if(myid % 2 != 0
		) {
			lindex=0;
			for(i=myst; i<myend; i++) {
				vec[i] = result[lindex];
				lindex++;
			}
		}

	//	if(myid == 1) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 3) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}
	//	
	//	else if(myid == 5) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}
	//	
	//	else if(myid == 7) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 9) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 11) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 13) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 15) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 17) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 19) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 21) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 23) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 25) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 27) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 29) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}

	//	else if(myid == 31) {
	//		lindex=0;
	//		for(i=myst; i<myend; i++) {
	//			vec[i] = result[lindex];
	//			lindex++;
	//		}
	//	}


		// Barrier to let vec modify.
		pthread_barrier_wait (&barrier);
	}

	if(myid == 0) {
		cout << k << "\n";

		gettimeofday(&end, NULL);
		printf("time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec)
		  - (start.tv_sec * 1000000 + start.tv_usec)));
		for(i=0; i<vr-1; i++)
			cout << "\n" << i << " " << vec[i];
	}
}

void mat_normal() {
	double sum=0;
	for(int i=0; i<vr-1; i++) {
		sum = 0;
		for(int j=row_ptrs[i]; j<row_ptrs[i+1]; j++) {
			sum += vals[j];
		}
		for(int j=row_ptrs[i]; j<row_ptrs[i+1]; j++) {
			vals[j] = vals[j]/sum;
		}
	}
}

void mat_transpose(std::vector<int> &nrows, std::vector<int> &ncols, std::vector<double> &nvals) {
	int i,j,k,l;

	for(i=0; i<row_ptrs.size(); i++) 
		ncols.push_back(0);

	for (i=0; i<col_inds.size(); i++)
		ncols[col_inds[i] + 1]++;
	for (i=0; i<ncols.size()-1; i++)
		ncols[i+1] += ncols[i];

	int *ptr = &row_ptrs.front();
	for (i=0; i<row_ptrs.size()-1; i++) {
		for(j=*ptr; j<*(ptr+1); j++) {
			k = col_inds[j];
			l = ncols[k]++;
			nrows[l] = i;
			nvals[l] = vals[j];
		}
		ptr++;
	}
	for (i = ncols.size()-1; i>0; i--)
		ncols[i] = ncols[i - 1];
	ncols[0] = 0;
}


int main(int argc, char *argv[]) {
	int i,j;

	// Reading graphs.
	read_graphs(argv);

	// Normalizing the matrix.
	mat_normal();

	// Transpose of the matrix.
	std::vector<int> nrows, ncols;
	std::vector<double> nvals;
	ncols.resize(row_ptrs.size());
	nrows.resize(col_inds.size());
	nvals.resize(vals.size());

	mat_transpose(nrows, ncols, nvals);

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

	for(i=0; i<vr; i++) 
		row_ptrs[i] = ncols[i];
	for(i=0; i<size; i++) 
		col_inds[i] = nrows[i];
	for(i=0; i<size; i++)	
		vals[i] = nvals[i];

	// Creating partitions.	
	for(i=0; i<partitions; i++) 
		part[i] = 0;
	j = (vr-1)%partitions;
	for(i=1; i<partitions; i++) {
		part[i] = part[i-1] + (vr-1)/partitions;
	}

//	int numberOfProcessors = sysconf(_SC_NPROCESSORS_ONLN);
//	cout << "Number of processors: " << numberOfProcessors;

	pthread_t p_threads[partitions]; 
	pthread_barrier_init(&barrier, NULL, partitions);	
	pthread_mutex_init(&minimum_value_lock, NULL);

//	pthread_attr_t attr;
//	cpu_set_t cpus;
//	pthread_attr_init(&attr);

	for (i=0; i<partitions; i++) {
		int *arg_send;
		arg_send = new int(3);
		arg_send[0] = i;			// Thread id.
		arg_send[1] = part[i];			// Thread index start.
		if(i != partitions-1 ) {		// Thread index end.
			arg_send[2] = part[i+1];
		}
		else {
			arg_send[2] = part[i]+((vr-1)/partitions)+j;	// All the remaining elements.
		}

//		CPU_ZERO(&cpus);
//		CPU_SET(i, &cpus);
//		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);
//		pthread_create(&p_threads[i], &attr, perform_matvec, (void *)arg_send);
		pthread_create(&p_threads[i], NULL, perform_matvec, (void *)arg_send);
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
