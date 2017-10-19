mpicc -O -o hmd_1thrd hmd_1.c -lm -fopenmp
mpirun --bind-to none -np 1 -machinefile $PBS_NODEFILE ./hmd_1thrd

mpicc -O -o hmd_2thrd hmd_2.c -lm -fopenmp
mpirun --bind-to none -np 1 -machinefile $PBS_NODEFILE ./hmd_2thrd

mpicc -O -o hmd_4thrd hmd_3.c -lm -fopenmp
mpirun --bind-to none -np 1 -machinefile $PBS_NODEFILE ./hmd_4thrd

mpicc -O -o hmd_8thrd hmd_4.c -lm -fopenmp
mpirun --bind-to none -np 1 -machinefile $PBS_NODEFILE ./hmd_8thrd
