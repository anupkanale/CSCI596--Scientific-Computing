#include "mpi.h"
#include <stdio.h>
#define NBIN 100000000

int nprocs;  /* Number of processors */
int myid;    /* My rank */

double global_sum(double partial) {
	double mydone, hisdone, partner;
	int bitValue;
	MPI_Status status;
	
  mydone = partial;
  for(bitValue=1; bitValue<nprocs; bitValue*=2){
  	partner=myid^bitValue;
 		MPI_Send(&mydone,1,MPI_DOUBLE,partner,bitValue,MPI_COMM_WORLD);
		MPI_Recv(&hisdone,1,MPI_DOUBLE,partner,bitValue,MPI_COMM_WORLD, &status);
    mydone = mydone + hisdone;
  }
  return mydone;
}

int main(int argc, char *argv[]) {
	int i;
	double step, x, partialSum=0.0, pi, sum=0.0;
	double cpu1, cpu2, elTime;
	step = 1.0/NBIN;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	cpu1 = MPI_Wtime();
	for (i=myid; i<NBIN; i+=nprocs) {
		x = (i+0.5)*step;
		partialSum += 4.0/(1.0+x*x);
	}

	partialSum *= step;
	pi = global_sum(partialSum);

/*
	for (i=0; i<NBIN; i++) {
		x = (i+0.5)*step;
		sum += 4.0/(1.0+x*x);
	}
	pi = sum*step;
*/

	cpu2 = MPI_Wtime();
	elTime = cpu2-cpu1;

	if (myid==0){
		printf("PI = %le\n", pi);
		printf("Time = %le\n", elTime);
	}

	MPI_Finalize();
	return 0;
}
