#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
 int main(int argc, char *argv[]) {
    // set up MPI
    MPI_Init(&argc, &argv);
    // get communicator size and my rank
    MPI_Comm comm = MPI_COMM_WORLD;
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);
    /* code */
    int N, C;
    // send N, C to the master processor 0
    if (rank == 0)
    {
        N = atoi(argv[1]);
        C = atoi(argv[2]);
    }

    // broadcast N, C to all the processors
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&C, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double time;
    time = MPI_Wtime(); // time start
    double sum = 0;

     // generate random numbers and calculate local sum

    double A[N/p]; // to store the random numbers
    sum = 0;
    srand48(C+rank); // set random seed
    for (int i=0; i < N/p; i++)
    {
        A[i] = drand48(); // generate random number
    }

    for (int i=0; i < N/p; i++)
    {
        sum += A[i]; // compute local sum
    }
    
    // find log(p)
    int d = 0;
    while((1<<d) < p)
    {
        d += 1;
    }
    
    // send and recieve messages
    for (int j=0; j < d; j++)
    {
        if ((rank & 1<<j) != 0)
        {
            MPI_Send(&sum, 1, MPI_DOUBLE, rank ^ (1<<j), 111, MPI_COMM_WORLD);
            return MPI_Finalize();
        }
        else 
        {
            double sum2;
            MPI_Status stat;
            MPI_Recv(&sum2, 1, MPI_DOUBLE, rank ^ (1<<j), 111, MPI_COMM_WORLD, &stat);

            sum += sum2; // add the received sum
        }

    }   
    
   
    time = MPI_Wtime() - time; // calculate runtime

    // print the results in the master processor 0
    if (rank == 0)
    {
    	printf("N = %i, P = %i, C = %i, S = %f \n", N, p, C, sum); 
    	printf("Time = %f  \n", time);

    }

    // finalize MPI
    MPI_Finalize();
    return 0; }





      
