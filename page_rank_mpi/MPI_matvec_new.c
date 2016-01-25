#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
#include<string.h>
#include <mpi.h>

void read_graphs(int argc, char *argv[], int *size, int *vr, int *vals, int *col_inds, int *row_ptrs, int *vec, int partitions, int *store_partitions[16], int *part) {
	//printf("yoyoyo");
	long i,j;
	int x; 
	ssize_t bytes_read;
	size_t nbytes=1000;
	char *mybuffer = NULL;

	// Reading the graph file.
	FILE *fp;
	//printf("%s ",argv[1]);
	fp = fopen(argv[1],"r");

	if(fp == NULL)
		exit(EXIT_FAILURE);
	
	i=0;
	char *tken;
	if((bytes_read = getline (&mybuffer, &nbytes, fp)) != -1) {
		tken = strtok(mybuffer," ");
		while(tken != NULL) {
			//printf("%s \n ",tken);
			vals[i] = atoi(tken);
			tken = strtok(NULL," ");
			i++;
		}
	}
	*size = i;

	i=0;
	if((bytes_read = getline (&mybuffer, &nbytes, fp)) != -1) {
		tken = strtok(mybuffer," ");
		while(tken != NULL) {
			//printf("%s \n ",tken);
			col_inds[i] = atoi(tken);
			tken = strtok(NULL," ");
			i++;
		}
	}

	i=0;
	if((bytes_read = getline (&mybuffer, &nbytes, fp)) != -1) {
		tken = strtok(mybuffer," ");
		while(tken != NULL) {
			//printf("%s \n ",tken);
			row_ptrs[i] = atoi(tken);
			tken = strtok(NULL," ");
			i++;
		}
	}
	*vr = i;

	fclose(fp);

	for(i=0; i<partitions; i++) {
		part[i] = 0;
	}


	fp = fopen(argv[2],"r");						// Time to read from the partition file
	if(fp == NULL)
		exit(EXIT_FAILURE);
	
	i=0;
	while((bytes_read = getline (&mybuffer, &nbytes, fp)) != -1) {
		x = atoi(mybuffer);
		store_partitions[x][part[x]] = i;
		part[x] = part[x]+1;
		i++;
	}

	//printf("%ld \n",i);

	fclose(fp);


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
//
//	for(i=0; i<partitions; i++) {
//		printf("\n Partitions %ld:  ",i);
//		for(j=0; j<part[i]; j++) {
//			printf("%d ",store_partitions[i][j]);
//		}
//	}

//	fflush(stdout);

	//Generating Vector.	//printf("Enter vector: \n");
	for(i=0; i<*vr-1; i++) {
		vec[i] = 1;
	}
}


void matvec(int myrank, long totElems, long nvr, long *mvals, long *mcols, long *mrows, double *mvecs, 
		long *rvmap, double *result, int partitions, int *store_partitions[16], int *local_partitions, int *part) {
	long i,j,k,index,kk,l=0;
	int flag;

	//printf("\n reached matvec");

	long *columns_needed = malloc(totElems * sizeof(long*));
	int total_columns_needed=0;
	for(j=0; j<totElems; j++) {
		flag = 0;
		for(kk=0; kk<nvr; kk++) {
			if(mcols[j] == rvmap[kk]) {
				flag = 1;
				break;
			}
		}
		if(flag == 0) {
			for(i=0; i<total_columns_needed; i++) {
				if(mcols[j] == columns_needed[i]) {
					flag = 1;
					break;
				}
			}
			if(flag == 0) {
				columns_needed[total_columns_needed] = mcols[j];
				total_columns_needed++;
			}
		}
	}

//	printf("\n Partions %d !! \n",myrank);
//	for(i=0; i<total_columns_needed; i++) {
//		printf("%ld -- ", columns_needed[i]);
//	}
//	printf("\n");

	// Now we need to find which partitions consists of which columns.
	long *part_columns[16];
	int *npart 	= malloc(partitions * sizeof(int));
	int *srdispls 	= malloc(partitions *sizeof(int));
	int *sdispls 	= malloc(partitions *sizeof(int));
	int *rdispls 	= malloc(partitions *sizeof(int));
	int *recvbuf 	= malloc(partitions *sizeof(int));	// This tells how many elements each partition needs.
	int *srcounts 	= malloc(partitions *sizeof(int));
	
	for(i=0; i<partitions; i++) {
		npart[i] = 0;
		part_columns[i] = malloc(total_columns_needed*sizeof(int));
	}

//	kk=0; flag = 0;								// For tracking elements to send.
//	for(j=0; j<partitions; j++) {						// For all the partitions.
//		if(j != myrank) {
//			flag = 0;
//			for(i=0; i<total_columns_needed; i++) {			// For all the columns needed.
//				for(k=0; k<part[j]; k++) {			// For all the enteries.
//					if(store_partitions[j][k] == columns_needed[i]) {
//						//printf("\n yes...");
//						part_columns[j][npart[j]] = columns_needed[i];
//						npart[j] = npart[j]+1;
//						flag = 1;	kk++;
//						break;
//					}
//				}
//			}
//			if(flag == 0) {
//				//printf("\n Myrank: %d -- j: %ld",myrank,j);
//				kk++;
//			}
//		}
//		else {
//			kk++;
//		}
//	}

	j=0; kk=0; flag=0;
	for(i=0; i<total_columns_needed; i++) {
		flag = 0;
		for(j=0; j<partitions; j++) {
			if(j != myrank) {
				for(k=0; k<part[j]; k++) {
					if(store_partitions[j][k] == columns_needed[i]) {
						part_columns[j][npart[j]] = columns_needed[i];
						npart[j] = npart[j]+1;
						flag = 1;	kk++;
						break;
					}
				}
			}
			if(flag == 1) {
				break;
			}
		}
	}

	for(i=0; i<partitions; i++) {
		if(npart[i] == 0) {
			kk++;
		}
	}

	//fflush(stdout);

//	printf("\n Myrank :: %d -- Total: %d -- kk: %ld    ",myrank,total_columns_needed,kk);
//
//	// Lets get things printed..
//	for(j=0; j<partitions; j++) {
//		if(j != myrank) {
//			//printf("\n My rank: %d --> Partition: %ld  -- Part[j]: %ld  !! ",myrank,j,part[j]);
//			k = npart[j];
//			for(i=0; i<k; i++) {
//				printf(" %ld -- ",part_columns[j][i]);
//			}
//		}
//	}
//	fflush(stdout);

	/* So we will have communication in three steps for the first round. 
	 * First each partition will tell all the other partitions the number of nodes they 
	 * expect to be sent. 
	 * Next, they will send the column ids they want from each partition.
	 * Final Communication will be sending the column values.
	 * Except for the first iteration, all other iterations will have only final communication.
	 */
	
	//int *sendcounts = malloc(partitions *sizeof(int));
	for(i=0; i<partitions; i++) {
		srdispls[i] 	= i;
		srcounts[i] 	= 1;
	}

	// First communication.
	MPI_Alltoallv(npart, srcounts, srdispls, MPI_INT, recvbuf, srcounts, srdispls, MPI_INT, MPI_COMM_WORLD);

//	printf("\n");
//	for(i=0; i<partitions; i++) {
//		printf("Elements %ld want from me(%d) -- %d \n",i,myrank,recvbuf[i]);
//	}
//	fflush(stdout);

	total_columns_needed = kk;				// Updated value.
	kk = 0;
	long *sendbuf = malloc(total_columns_needed * sizeof(long));
	for(i=0; i<partitions; i++) {				// Now we are trying to place the columns at the locations.
		sdispls[i] = kk;
		if(i != myrank) {
			k = npart[i];
			for(j=0; j<k; j++) {
				sendbuf[kk] = part_columns[i][j];
				kk++;
			}
			if(j == 0) {				// If no element required from this partition.
				sendbuf[kk] = -1;
				kk++;
			}
		}
		else {
			sendbuf[kk] 	= -1;			// Receive from self a negative.
			kk++;
		}
	}

	//printf("\n value of total: %d -- kk: %ld ",total_columns_needed,kk);

//	printf("\n Columns needed by me(%d) -- ",myrank);
//	for(i=0; i<total_columns_needed; i++) {
//		printf("%ld ",sendbuf[i]);
//	}
//
//
//	printf("\n Displacement needed by me(%d) -- ",myrank);
//	for(i=0; i<partitions; i++) {
//		printf("%d ",sdispls[i]);
//	}


	int total_columns_send = 0;
	for(i=0; i<partitions; i++) {
		if(recvbuf[i] > 0) {
			total_columns_send += recvbuf[i];
		}
		else {						// For self or empty partitions
			total_columns_send++;
			recvbuf[i] = 1;				// TODO -- Basically for sending again.
		}
	}

	kk = 0;
	for(i=0; i<partitions; i++) {
		rdispls[i] = kk;
		if(i != myrank) {
			kk += recvbuf[i];
		}	
		else {
			kk++;
		}
	}

	long *receivebuf = malloc(total_columns_send *sizeof(long)); 

	// TODO -- We will now be modifying the part array by setting 1 for all 0, i.e. 
	// allowing transfer of on element.
	for(i=0; i<partitions; i++) {
		if(npart[i] == 0)
			npart[i] = 1;
	}

//	int temp=0;
//	for(i=0; i<partitions; i++) {
//		temp += npart[i];
//	}

	//int *rcounts = malloc(partitions *sizeof(int));
	//printf("\n I am partition %d -- My sendbuf is: %d (%d) -- My recvbuf is: %d \n",myrank,total_columns_needed,temp,total_columns_send);

      	// Second communication.
	MPI_Alltoallv(sendbuf, npart, sdispls, MPI_LONG, receivebuf, recvbuf, rdispls, MPI_LONG, MPI_COMM_WORLD);
	
//	j=0;
//	for(i=0; i<total_columns_send; i++) {
//		if(i == rdispls[j]) {
//			printf("\n I %d need to send following columns to partition %ld: ",myrank,j);
//			j++;
//		}
//		printf("%ld  ",receivebuf[i]);
//	}


	/* Now everyone knows what all columns of the vector they have to send to the requesting partitions.
	 * So everyone starts by collecting all the column enteries from their respective vectors and sending the same. 
	 * The point to note here is that vec (vector) is represented as a pointer. 
	 * Again we will be having an all-to-all communication, but now this time there will be a 
	 * converging loop, which will determine the extent of algorithm.
	 */

	double *sbuffer = malloc(total_columns_send * sizeof(double));		// Used for other calculation.
	double *rbuffer = malloc(total_columns_needed * sizeof(double));	// Used for my calculation.
	int *locmap 	= malloc(total_columns_send * sizeof(int));
	
	for(i=0; i<total_columns_send; i++) {
		sbuffer[i] = receivebuf[i];
		locmap[i] = 0;
	}

	for(i=0; i<nvr; i++) {
		result[i] = 0;
	}

	

	int converge = 1; 
	double negnum = -1, sum=0, temp=0, total_sum=0, check=0.00001;
	for(i=0; i<nvr; i++) {							// Computing the local vector sum.
		sum += mvecs[i];
	}

	flag = 0;
	for(i=0; i<total_columns_send; i++) {			// The list of values to send.
		flag = 0;
		if(receivebuf[i] != negnum) {
			k = receivebuf[i];
			for(j=0; j<nvr; j++) {
				if(rvmap[j] == k) {
					locmap[i] = j;
					flag = 1;
					break;
				}
			}
		}
		else {
			locmap[i] = -1;
			flag = 1;
		}
		if(flag == 0) {
			locmap[i] = -1;
		}
	}

	int xsum=0;
	for(i=0; i<nvr; i++) {
		for(j=mrows[i]; j<mrows[i+1]; j++) {
			xsum++;	
		}
	}
	int *locpos = malloc(xsum * sizeof(int));
	int *locbool = malloc(xsum * sizeof(int));
	
	l=0;
	for(i=0; i<nvr; i++) {
		for(j=mrows[i]; j<mrows[i+1]; j++) {
			k = mcols[j];					// Need to fetch this col.
			flag = 0;
			for(kk=0; kk<nvr; kk++) {
				if(rvmap[kk] == k) {
					locpos[l] = kk;
					locbool[l] = 0;
					l++;
					break;
				}
			}
			if(flag == 0) {
				for(kk=0; kk<total_columns_needed; kk++) {
					if(sendbuf[kk] == k) {
						locpos[l] = kk;
						locbool[l] = 1;
						l++;
						break;
					}
				}
			}
		}
	}


//	fflush(stdout);

	MPI_Allreduce(&sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	// Sum of all the vec elements.

	for(i=0; i<nvr; i++) {
		mvecs[i] = mvecs[i]/total_sum;
	}

	// Getting the start time.
	struct timeval start, end;
	if(myrank == 0) {
		gettimeofday(&start, NULL);
	}

	int ll=0;
	while(converge == 1) {
		sum = 0;

	//	for(i=0; i<total_columns_send; i++) {			// The list of values to send.
	//		if(receivebuf[i] != negnum) {
	//			k = receivebuf[i];
	//			for(j=0; j<nvr; j++) {
	//				if(rvmap[j] == k) {
	//					sbuffer[i] = mvecs[j];
	//					break;
	//				}
	//			}
	//		}
	//	}
	
		for(i=0; i<total_columns_send; i++) {
			if(locmap[i] != -1) {
				sbuffer[i] = mvecs[locmap[i]];
			}
		}

		// Third communication -- reverse of receivebuf and sendbuf.
		MPI_Alltoallv(sbuffer, recvbuf, rdispls, MPI_DOUBLE, rbuffer, npart, sdispls, MPI_DOUBLE, MPI_COMM_WORLD);

		// The matrix multiplication part.
		double valneed;
		l=0;
		for(i=0; i<nvr; i++) {
			for(j=mrows[i]; j<mrows[i+1]; j++) {
			//	k = mcols[j];					// Need to fetch this col.
			//	flag = 0;
			//	for(kk=0; kk<nvr; kk++) {
			//		if(rvmap[kk] == k) {
			//			valneed = mvecs[kk];
			//			flag = 1;
			//			break;
			//		}
			//	}
			//	if(flag == 0) {
			//		for(kk=0; kk<total_columns_needed; kk++) {
			//			if(sendbuf[kk] == k) {
			//				valneed = rbuffer[kk];
			//				break;
			//			}
			//		}
			//	}
			      
			//	result[i] += mvals[j] * valneed;
				if(locbool[l] == 1) {
					result[i] += mvals[j] * rbuffer[locpos[l]];
				}
				else {
					result[i] += mvals[j] * mvecs[locpos[l]];
				}
				l++;
			}
		}

		for(i=0; i<nvr; i++) {							// Local result summation.
			sum += result[i];
		}

		MPI_Allreduce(&sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	// Sum of all the vec elements.

		for(i=0; i<nvr; i++) {						
			result[i] = result[i]/total_sum;
		}
	
		sum = 0;
		for(i=0; i<nvr; i++) {
			temp = result[i] - mvecs[i];
			temp = pow(temp,2);
			sum += temp;
		}

		MPI_Allreduce(&sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	// Sum of all the vec elements.

		total_sum = sqrt(total_sum);
		if(total_sum < check) {
			converge = 0;
		}
		
		// Trasferring result to mvecs.
		for(i=0; i<nvr; i++) {
			mvecs[i] = result[i];
		}
	}

	// Getting the end time.
	if(myrank == 0) {
		gettimeofday(&end, NULL);
		printf("time: %ld\n", ((end.tv_sec * 1000000 + end.tv_usec)
			  - (start.tv_sec * 1000000 + start.tv_usec)));
	}
}


void mpi_pagerank(int myrank, int size, int vr, int *vals, int *col_inds, int *row_ptrs, int *vec, int partitions, int *store_partitions[16], int *part) {
	int *local_partitions; //= malloc(enteries*sizeof(int));

	// The appropriate size data structures.
	long totElems=0, nvr, temp;	
	long *mvals, *mcols, *mrows, *rvmap; 
	double *result, *mvecs;
	
	if(myrank == 0) {
		int i,j,k,kk;
		int *tg_partitions; 

		for(i=1; i<partitions; i++) {
			MPI_Send(&size,		1, 		MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&vr, 		1, 		MPI_INT, i, 1, MPI_COMM_WORLD);
			MPI_Send(&partitions, 	1, 		MPI_INT, i, 2, MPI_COMM_WORLD);
			MPI_Send(part, 		partitions, 	MPI_INT, i, 3, MPI_COMM_WORLD);

			nvr 	= part[i];
			mrows 	= malloc((nvr+1) * sizeof(long*));
        	        rvmap 	= malloc(nvr * sizeof(long*));
        	        mvecs 	= malloc(nvr * sizeof(double*));
			result 	= malloc(nvr * sizeof(double*));

			// Constructing the vector.
			for(j=0; j<nvr; j++) {
				mvecs[j] = 1;
			}

			// Sending all the partitions to each node.
			for(k=0; k<partitions; k++) {
				local_partitions = malloc(part[k]*sizeof(int));
				for(j=0; j<part[k]; j++) {
					local_partitions[j] = store_partitions[k][j];
				}
				if(k == i) {
					tg_partitions = malloc(part[k]*sizeof(int));
					for(j=0; j<part[k]; j++) {
						tg_partitions[j] = store_partitions[k][j];
					}
				}
				MPI_Send(local_partitions, part[k], MPI_INT, i, 4*partitions+k, MPI_COMM_WORLD);
			}

			// Creating the individual row pointer.
			long numElem = 0; kk=0; temp=0; totElems=0;
			for(j=0; j<part[i]; j++) {
				if(j != 0) {
					mrows[j] = numElem + mrows[j-1];		// Storing the row at its new position.
				}
				else {
					mrows[j] = numElem;
				}
				temp = tg_partitions[j]; 
				rvmap[j] = temp;
				numElem = row_ptrs[temp+1]-row_ptrs[temp];		// Getting number of elements in a row.
				for(k=row_ptrs[temp]; k<row_ptrs[temp+1]; k++) {
					totElems++;
				}
			}
			//printf("\n Total elems: %ld",totElems);
			mrows[j] = numElem + mrows[j-1];				// Last entry of new csr.

			mvals = malloc(totElems * sizeof(long*));
                	mcols = malloc(totElems * sizeof(long*));
			for(j=0; j<part[i]; j++) {
				temp = tg_partitions[j]; 
				for(k=row_ptrs[temp]; k<row_ptrs[temp+1]; k++) {				
					mvals[kk] = vals[k];				// Value subset created.
					mcols[kk] = col_inds[k];			// Column subset created.
					kk++;
				}
			}

		//	printf("\n For partition %d -- info is: ",i);
		//	printf("\n vals: ");
		//	for(k=0; k<totElems; k++) {
		//		printf("%ld ",mvals[k]); 
		//	}
		//	printf("\n cols: ");
		//	for(k=0; k<totElems; k++) {
		//		printf("%ld ",mcols[k]); 
		//	}
		//	printf("\n rows: ");
		//	for(k=0; k<nvr+1; k++) {
		//		printf("%ld ",mrows[k]); 
		//	}
		//	printf("\n rvmap: ");
		//	for(k=0; k<nvr; k++) {
		//		printf("%ld ",rvmap[k]); 
		//	}
		//	printf("\n");
		//	fflush(stdout);

			MPI_Send(&totElems, 	1, 		MPI_LONG, 	i, 5, MPI_COMM_WORLD);
			MPI_Send(mvals, 	totElems, 	MPI_LONG, 	i, 6, MPI_COMM_WORLD);
			MPI_Send(mcols, 	totElems, 	MPI_LONG, 	i, 7, MPI_COMM_WORLD);
			MPI_Send(mrows, 	nvr+1, 		MPI_LONG, 	i, 8, MPI_COMM_WORLD);
			MPI_Send(mvecs, 	nvr, 		MPI_DOUBLE, 	i, 9, MPI_COMM_WORLD);
			MPI_Send(rvmap, 	nvr, 		MPI_LONG, 	i, 10, MPI_COMM_WORLD);
		}

		// Setting the similar data structures for partition zero.
		local_partitions = malloc(part[0]*sizeof(int));
		for(j=0; j<part[0]; j++) {
			local_partitions[j] = store_partitions[myrank][j];
		}

		long numElem = 0; kk=0; temp=0; totElems=0;
		for(j=0; j<part[0]; j++) {
			if(j != 0) {
				mrows[j] = numElem + mrows[j-1];		// Storing the row at its new position.
			}
			else {
				mrows[j] = numElem;
			}
			temp = local_partitions[j]; 
			rvmap[j] = temp;
			numElem = row_ptrs[temp+1]-row_ptrs[temp];		// Getting number of elements in a row.
			for(k=row_ptrs[temp]; k<row_ptrs[temp+1]; k++) {
				totElems++;
			}
		}
		mrows[j] = numElem + mrows[j-1];				// Last entry of new csr.

		mvals = malloc(totElems * sizeof(long*));
                mcols = malloc(totElems * sizeof(long*));
		for(j=0; j<part[0]; j++) {
			temp = local_partitions[j]; 
			for(k=row_ptrs[temp]; k<row_ptrs[temp+1]; k++) {				
				mvals[kk] = vals[k];				// Value subset created.
				mcols[kk] = col_inds[k];			// Column subset created.
				kk++;
			}
		}

		free(vals);	free(col_inds);		free(row_ptrs);		free(vec);
	}
	else {
		int i,j,k;
		//printf("My rank: %d\n",myrank);
		
		MPI_Recv(&size, 	1, 		MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&vr, 		1, 		MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&partitions, 	1, 		MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		part = malloc(partitions * sizeof(int*));
		MPI_Recv(part, 		partitions, 	MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		nvr 	= part[myrank];
		mrows 	= malloc((nvr+1) * sizeof(long*));
		rvmap 	= malloc(nvr * sizeof(long*));
		mvecs 	= malloc(nvr * sizeof(double*));
		result 	= malloc(nvr * sizeof(double*));

		// Initializing the all partition array.
		for(i=0; i<partitions; i++) {
			store_partitions[i] = malloc(part[i]*sizeof(int));
		}

		// Receiving all the partitions.
		for(k=0; k<partitions; k++) {
			local_partitions = malloc(part[k]*sizeof(int));
			MPI_Recv(local_partitions, part[k], MPI_INT, 0, 4*partitions+k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
		//	printf("My rank: %d\n",myrank);
		//	printf("\n Printing Local Partition \n");
		//	for(j=0; j<part[k]; j++) {
		//		printf("%d ",local_partitions[j]);
		//	}
		//	fflush(stdout);

			for(j=0; j<part[k]; j++) {
				store_partitions[k][j] = local_partitions[j];
			}
		}

		// Storing the local partition.
		local_partitions = malloc(part[myrank]*sizeof(int));
		for(j=0; j<part[myrank]; j++) {
			local_partitions[j] = store_partitions[myrank][j];
		}
		
	//	printf("My rank: %d\n",myrank);
	//	printf("\n Printing Local Partition \n");
	//	for(i=0; i<part[myrank]; i++) {
	//		printf("%d ",local_partitions[i]);
	//	}
	//	fflush(stdout);

		totElems = 0;
		MPI_Recv(&totElems, 1, MPI_LONG, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		//printf("\n My rank:  %d -- Tot elems: %ld",myrank,totElems);

		mvals = malloc(totElems * sizeof(long*));
		mcols = malloc(totElems * sizeof(long*));

		//printf("\n My rank: %d -- nvr: %ld  -- %d",myrank,nvr+1,vr);

		MPI_Recv(mvals, totElems, 	MPI_LONG, 	0, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(mcols, totElems, 	MPI_LONG, 	0, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(mrows, nvr+1, 		MPI_LONG, 	0, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(mvecs, nvr, 		MPI_DOUBLE, 	0, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(rvmap, nvr, 		MPI_LONG, 	0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		long ii;
	//	printf("My rank: %d\n",myrank);
	//	printf("\n Vals: \n");
	//	for(ii=0; ii<totElems; ii++) {
	//		printf("%ld ",mvals[ii]);
	//	}
	//
	//	printf("\n Cols: \n");
	//	for(i=0; i<totElems; i++) {
	//		printf("%ld ",mcols[i]);
	//	}
	//	printf("\n");
	//
	//	printf("\n Rows: \n");
	//	for(i=0; i<nvr+1; i++) {
	//		printf("%ld ",mrows[i]);
	//	}
	//	printf("\n");

	//	printf("\n Rows Map: \n");
	//	for(i=0; i<nvr; i++) {
	//		printf("%ld ",rvmap[i]);
	//	}
	//	printf("\n");
	//	fflush(stdout);

	//	printf("\n Done \n");
	}


	MPI_Barrier(MPI_COMM_WORLD);	

	if(myrank == 0) {
		matvec(myrank,totElems,nvr,mvals,mcols,mrows,mvecs,rvmap,result,partitions,store_partitions,local_partitions,part);
	}
	else {
		matvec(myrank,totElems,nvr,mvals,mcols,mrows,mvecs,rvmap,result,partitions,store_partitions,local_partitions,part);
	}

	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	// Everyone prints the result.
	if(myrank == 0) {
		printf("node_Id pagerank \n");
	}
	else {
		;
	}

	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	int i=0,j=0;
//	while(j<vr) {
//		if(j == rvmap[i]) {
//			printf("%ld %f\n",rvmap[i],mvecs[i]);
//			i++;
//		}
//		j++;
//		fflush(stdout);
//		MPI_Barrier(MPI_COMM_WORLD);
//	}

	if(myrank == 0) {
		for(i=0; i<nvr; i++) {
			printf("%ld %0.23f\n",rvmap[i],mvecs[i]);
		}
		fflush(stdout);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(myrank == 1) {
		//printf("!!!!!!!!!!!!!!");
		for(i=0; i<nvr; i++) {
			printf("%ld %0.23f\n",rvmap[i],mvecs[i]);
		}
		fflush(stdout);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(myrank == 2) {
		//printf("!!!!!!!!!!!!!!");
		for(i=0; i<nvr; i++) {
			printf("%ld %0.23f\n",rvmap[i],mvecs[i]);
		}
		fflush(stdout);	
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	if(myrank == 3) {
		//printf("!!!!!!!!!!!!!!");
		for(i=0; i<nvr; i++) {
			printf("%ld %0.23f\n",rvmap[i],mvecs[i]);
		}
		fflush(stdout);
	}

//	for(i=0; i<nvr; i++) {
//		printf("%ld %f\n",rvmap[i],mvecs[i]);
//	}
}



int main(int argc, char *argv[]) {
	int numprocs, myrank;

	int size, vr, partitions;
	int *vals 	;//= malloc(4000000*sizeof(int));
	int *col_inds 	;//= malloc(4000000*sizeof(int));
	int *row_ptrs 	;//= malloc(4000000*sizeof(int));
	int *vec 	;//= malloc(4000000*sizeof(int));
	int *part;
	//int *result 	= malloc(4000000*sizeof(int));
	int *store_partitions[16]; 

	

	MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 

	if(myrank == 0) {
		partitions = atoi(argv[3]);							// Command line argument specifies the number of partitions.
		part = malloc(partitions * sizeof(int*));

		vals 		= malloc(80000000*sizeof(int));
                col_inds 	= malloc(80000000*sizeof(int));
                row_ptrs 	= malloc(80000000*sizeof(int));
                vec 		= malloc(80000000*sizeof(int));

		int ii;
		for(ii=0; ii<partitions; ii++) {
			store_partitions[ii] = malloc(60000000*sizeof(int));
		}

		read_graphs(argc,argv,&size,&vr,vals,col_inds,row_ptrs,vec,partitions,store_partitions,part);
		
	}
	else {
		;	// Do nothing
	}

	MPI_Barrier(MPI_COMM_WORLD);		

	mpi_pagerank(myrank,size,vr,vals,col_inds,row_ptrs,vec,partitions,store_partitions,part);
	
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
}

