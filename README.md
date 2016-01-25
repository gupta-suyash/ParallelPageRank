# ParallelPageRank

This work dealt with developing the Parallel implementation of Page Rank Algorithm. We implemented the Page Rank algorithm using the Matrix-Vector multiplication and Random Walk. For the matrix-vector product technique we computed the L1 normalization to converge the page rank. We implemented this technique using both MPI and Pthreads. The random walk approach used the random number distribution present in Pthreads. The algorithms were able to provide us an efficiency of 66% for the 4 million nodes graph.
