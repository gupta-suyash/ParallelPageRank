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
#include<sys/time.h>

using namespace std;

// Defining global data.
const int partitions = 1;			// Determines the number of threads to run. 
int size, vr, part[partitions];
double total_sum[partitions];
std::vector<int> col_inds, row_ptrs;
std::vector<double> vec, vals;


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


void perform_matvec() {
	long i,j,flag=1, k=0, lindex;
	double temp, sum, check = 0.0000001;

	std::vector<double> result;
        result.resize(vr-1);

	for(i=0; i<vr-1; i++) 
		result[i] = 0;

	struct timeval start, end;
	gettimeofday(&start, NULL);

	//for(int ll=0; ll<1; ll++) {
	int iter=0;
	for(iter=0; iter<8; iter++)
	//while(flag == 1) 
	{
		sum = 0;	temp=0;
		lindex=0;
		for(i=0; i<vr-1; i++) {
			temp = 0;
			for(j=row_ptrs[i]; j<row_ptrs[i+1]; j++) 
				temp += vals[j] * vec[col_inds[j]];
			result[lindex] = temp;
			lindex++;
		}

		lindex=0;
		for(i=0; i<vr-1; i++) {
			temp = result[lindex] - vec[i];
			temp = pow(temp,2);
			sum += temp;
			lindex++;
		}

		sum = sqrt(sum);
		if(sum < check) {
			flag = 0;
		}

		k++;

	//	for(i=0; i<vr-1; i++)
	//		cout << "\n" <<  result[i];

		lindex = 0;
		for(i=0; i<vr-1; i++) {
			vec[i] = result[lindex];
			lindex++;
		}
	}

	//cout << k << "\n";

	gettimeofday(&end, NULL);
	printf("time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec)
	  - (start.tv_sec * 1000000 + start.tv_usec)));
	for(i=0; i<vr-1; i++)
		cout << "\n" << i << " " << vec[i];
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

	// Serial Matvec.
	perform_matvec();

	
//	// Printing Result.
//	printf("\n Result: \n");
//	for(i=0; i<vr-1; i++) {
//		printf("%d ",result[i]);
//	}
//
//	printf("\n Bye \n");
}
