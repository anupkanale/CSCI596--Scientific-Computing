#include "mpi.h"
#include <stdio.h>
int nprocs;  /* Number of processors */
int myid;    /* My rank */

double global_sum(double partial) {
  /* Implement your own global summation here */
	double hisdone, mydone, partner;
	int bitValue;
	MPI_Status status;
	
  mydone = partial;
  for (bitValue=1; bitValue<nprocs; bitValue*=2){
		//partner := myid XOR 2-to-the-power-l;
  	partner=myid^bitValue;
    //send mydone to partner;
		MPI_Send(&mydone,1,MPI_DOUBLE,partner,bitValue,MPI_COMM_WORLD);
    //receive hisdone from partner;
		MPI_Recv(&hisdone,1,MPI_DOUBLE,partner,bitValue,MPI_COMM_WORLD, &status);
    mydone = mydone + hisdone;
  }
  return mydone;
}

int main(int argc, char *argv[]) {
  double partial, sum, avg;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  partial = (double) myid;
  printf("Node %d has %le\n", myid, partial);

  sum = global_sum(partial);

  if (myid == 0) {
    avg = sum/nprocs;
    printf("Global average = %le\n", avg);
  }

  MPI_Finalize();

  return 0;
}
