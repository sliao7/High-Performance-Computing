#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string.h>
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

    srand48(C+rank); // set random seed

    double time;
    time = MPI_Wtime(); // start time
    double sum = 0;

     // generate random numbers and calculate local sum

    double A[N/p];
    sum = 0;
    
    for (int i=0; i < N/p; i++)
    {
        A[i] = drand48(); // generate random number
    }

    for (int i=0; i < N/p; i++)
    {
        sum += A[i]; // compute local sum
    }
    
    double local_time = MPI_Wtime() - time;
    time = MPI_Wtime();

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
            //return 0;
        }
        else 
        {
            double sum2;
            MPI_Status stat;
            MPI_Recv(&sum2, 1, MPI_DOUBLE, rank ^ (1<<j), 111, MPI_COMM_WORLD, &stat);

            sum += sum2;
        }

    }   
    
   
    double paralle_time = MPI_Wtime() - time; // running time


    if (rank == 0)
    {
    	printf("N = %i, P = %i, C = %i, S = %f \n", N, p, C, sum);
    	printf("%f, %f  \n", local_time, paralle_time);
        // printf("%f\n",time);
    }

    

   // printf("Hello from rank %i/%i.\n", rank, p);
    // finalize MPI
    MPI_Finalize();
    return 0; }





      
