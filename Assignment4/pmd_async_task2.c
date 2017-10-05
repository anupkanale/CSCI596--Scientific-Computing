/*----------------------------------------------------------------------
Program pmd.c performs parallel molecular-dynamics for Lennard-Jones 
systems using the Message Passing Interface (MPI) standard.
----------------------------------------------------------------------*/
#include "pmd_async_task2.h"
#define VMAX 5.0  // Max. velocity value to construct a velocity histogram
#define NBIN 100  // # of bins in the histogram

FILE *fpv;

/*--------------------------------------------------------------------*/
int main(int argc, char **argv) {
/*--------------------------------------------------------------------*/
  double cpu1,cpu;
	int i,a;

  MPI_Init(&argc,&argv); /* Initialize the MPI environment */
  //MPI_Comm_rank(MPI_COMM_WORLD, &sid);  /* My processor ID */

	MPI_Comm_rank(MPI_COMM_WORLD, &gid); //Global rank
	md_shit = gid%2; // = 1 (MD workers) or 0 (analysis workers)
	MPI_Comm_split(MPI_COMM_WORLD,md_shit,0,&workers);
	MPI_Comm_rank(workers,&sid); // Rank in workers

  /* Vector index of this processor */
  vid[0] = sid/(vproc[1]*vproc[2]);
  vid[1] = (sid/vproc[2])%vproc[1];
  vid[2] = sid%vproc[2];

	init_params();
	if (md_shit) {
		set_topology();
		init_conf();
		atom_copy();
		compute_accel();
	}
	else
		if (sid == 0)
		 fpv = fopen("pv.dat","w");

  //printf("%d\n", md_shit);
  cpu1 = MPI_Wtime();
  for (stepCount=1; stepCount<=StepLimit; stepCount++) {
	if (md_shit){
	 printf("%d\n", md_shit);
	 single_step();
	}
	if (stepCount%StepAvg == 0) {
		if (md_shit) {
			// Send # of atoms, n, to rank gid-1 in MPI_COMM_WORLD
			MPI_Send(&n, 1, MPI_INT, gid-1, 1000, MPI_COMM_WORLD);
			// Compose message to be sent 
			for(i=0; i<n; i++)
				for(a=0; a<3; a++)
					dbuf[3*i+a] = rv[i][a];
			// Send velocities of n atoms to rank gid-1 in MPI_COMM_WORLD
			MPI_Send(dbuf, 3*n, MPI_DOUBLE, gid-1, 2000, MPI_COMM_WORLD);
			eval_props();
		}
		else {
			// Receive # of atoms, n, from rank gid+1 in MPI_COMM_WORLD
			MPI_Recv(&n, 1, MPI_INT, gid+1, 1000, MPI_COMM_WORLD,&status);
			// Receive velocities of n atoms from rank gid+1 in MPI_COMM_WORLD
			MPI_Recv(dbufr, 3*n, MPI_DOUBLE, gid+1, 2000, MPI_COMM_WORLD,&status);
			for(i=0; i<n; i++)
				for(a=0; a<3; a++)
					rv[i][a] = dbufr[3*i+a];
			calc_pv();
		} // end if
	} // end if
  } // end for


  cpu = MPI_Wtime() - cpu1;
	if(md_shit && sid==0)
		printf("CPU & COMT = %le %le\n",cpu,comt);
	if(!md_shit && sid==0)
		fclose(fpv);

  MPI_Finalize(); /* Clean up the MPI environment */
  return 0;
}
