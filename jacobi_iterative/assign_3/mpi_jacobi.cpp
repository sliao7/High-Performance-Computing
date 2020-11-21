/**
 * @file    mpi_jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements MPI functions for distributing vectors and matrixes,
 *          parallel distributed matrix-vector multiplication and Jacobi's
 *          method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "mpi_jacobi.h"
#include "jacobi.h"
#include "utils.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <stdio.h>

/*
 * TODO: Implement your solutions here
 */


void distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
{
    int nlocal;
    int npes;
    int my2drank, my2dcoords[2];
    int coords[2];


    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &my2drank);
    MPI_Cart_coords(comm, my2drank, 2, my2dcoords);


    // get the number of elements that the current processor should receive
    int q = sqrt(npes);
    if (my2dcoords[1] == 0)
    {
        nlocal = block_decompose(n,q,my2dcoords[0]);
    }
    else
    {
        nlocal = 0;
    }

    *local_vector = new double[nlocal];

    int send_counts[npes];


    for (int i = 0; i < npes; i++)
    {
        MPI_Cart_coords(comm, i, 2, coords);
        if (coords[1] == 0)
        {
            send_counts[i] = block_decompose(n, q, coords[0]);
        }
        else
        {
            send_counts[i] = 0;
        }
    }

    int displ[npes];
    displ[0] = 0;
    for (int i = 1; i < npes; i++)
    {
        displ[i] = displ[i-1] + send_counts[i-1];
        
    }
    
    

    MPI_Scatterv(input_vector, send_counts, displ, MPI_DOUBLE, *local_vector, nlocal, MPI_DOUBLE, 0, comm);
    

 
}


// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm)
{
    // TODO
    int nlocal;
    int npes;
    int my2drank, my2dcoords[2];
    int coords[2];


    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &my2drank);
    MPI_Cart_coords(comm, my2drank, 2, my2dcoords);



    int q = sqrt(npes);

    // get the number of elements that the current processor should send
    if (my2dcoords[1] == 0)
    {
        nlocal = block_decompose(n,q,my2dcoords[0]);
    }
    else
    {
        nlocal = 0;
    }

    // get the receive counts
    int receiv_counts[npes];
    for (int i = 0; i < npes; i++)
    {
        MPI_Cart_coords(comm, i, 2, coords);
        if (coords[1] == 0)
        {
            receiv_counts[i] = block_decompose(n, q, coords[0]);
        }
        else
        {
            receiv_counts[i] = 0;
        }
    }

    int displ[npes];
    displ[0] = 0;
    for (int i = 1; i < npes; i++)
    {
        displ[i] = displ[i-1] + receiv_counts[i-1];
        
    }

    MPI_Gatherv(local_vector, nlocal, MPI_DOUBLE, output_vector, receiv_counts, displ, MPI_DOUBLE, 0, comm);

    
}


void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
{
    // 
    
    int nlocal;
    int npes;
    int my2drank, my2dcoords[2];
    int coords[2];


    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &my2drank);
    MPI_Cart_coords(comm, my2drank, 2, my2dcoords);

    int x = my2dcoords[0];
    int y = my2dcoords[1];

    int q = sqrt(npes);

    int nlocal_x = block_decompose(n,q,x);// number of rows of the lcoal_matrix

    //step 1: distribute the matrix elements into processors in the first column
    if (my2dcoords[1] == 0)
    {
        nlocal = n*nlocal_x;
    }
    else
    {
        nlocal = 0;
    }
    


    int send_counts[npes];


    for (int i = 0; i < npes; i++)
    {
        MPI_Cart_coords(comm, i, 2, coords);
        if (coords[1] == 0)
        {
            send_counts[i] = (block_decompose(n,q,coords[0]))*n;
        }
        else
        {
            send_counts[i] = 0;
        }
    }



    int displ[npes];
    displ[0] = 0;
    for (int i = 1; i < npes; i++)
    {
        displ[i] = displ[i-1] + send_counts[i-1];
        
    }
    
   
    double* local_vector = NULL;
    if(y == 0)
    {
        local_vector = new double[nlocal];
    }
    
   
    MPI_Scatterv(input_matrix, send_counts, displ, MPI_DOUBLE, local_vector, nlocal, MPI_DOUBLE, 0, comm);

    // Step2: distribute the numbers from processors in the first column to other processors in the same row
    // set up the subcommunications, processors in the same row goes to the same subcommunication
    int rank;
    MPI_Comm row_comm;
    MPI_Comm_split(comm, my2dcoords[0],my2drank,&row_comm);
    MPI_Comm_rank(row_comm, &rank);

    // arrange the local_vectors so that after scatterv we get the correct elements for local_matrices
    double* arranged = new double[nlocal];
    int index = 0;
    for(int j = 0; j < n; j++)
    {
        for(int i = 0; i < nlocal/n; i++)
        {
            arranged[index] = local_vector[n*i+j];
            index++;
        }
    }

    if (y == 0)
    {
        delete[] local_vector;
        local_vector = NULL;
    }
        

    int row_npes;
    MPI_Comm_size(row_comm, &row_npes);
    nlocal = (block_decompose(n,q,my2dcoords[0]))*n;
    int row_send_counts[row_npes];
    for (int i = 0; i < row_npes; i++)
    {      
        row_send_counts[i] = block_decompose(n,q,i)*nlocal_x;             
    }


    int row_displ[row_npes];
    row_displ[0] = 0;
    for (int i = 1; i < row_npes; i++)
    {
        row_displ[i] = row_displ[i-1] + row_send_counts[i-1];     
    }


    int row_nlocal = row_send_counts[rank]; // number of elements for each processors to receive 
    *local_matrix = new double[row_nlocal];
    double* local_matrix_arranged = new double[row_nlocal];

    MPI_Scatterv(arranged, row_send_counts, row_displ, MPI_DOUBLE, local_matrix_arranged, row_nlocal, MPI_DOUBLE, 0, row_comm);
    MPI_Comm_free(&row_comm);

    delete[] arranged;
    arranged = NULL;

    // arrange the received elements in the correct order of local matrices
    index = 0;
    int num_rows = nlocal/n;
    int num_columns = row_nlocal/(nlocal/n);

    for(int j = 0; j < num_rows; j++)
    {
        for(int i = 0; i < num_columns; i++)
        {
            (*local_matrix)[index] = local_matrix_arranged[num_rows*i+j];
            index++;
        }
    } 

    delete[] local_matrix_arranged;
    local_matrix_arranged = NULL;
     

}


void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
{
    // TODO
    int nlocal_x, nlocal_y;
    int npes;
    int my2drank, my2dcoords[2];


    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &my2drank);
    MPI_Cart_coords(comm, my2drank, 2, my2dcoords);
    MPI_Status status;

    int x = my2dcoords[0];
    int y = my2dcoords[1];

    int q = sqrt(npes);
   
    nlocal_x = block_decompose(n,q,x);
    nlocal_y = block_decompose(n,q,y);
    int reverse_coord[2] = {y,x};
    int reverse_rank;
    MPI_Cart_rank(comm, reverse_coord, &reverse_rank);

    // Step1: send the elements in processors with coord = [i,0] to processors with coord = [0,i],i.e., processors in the first column 
    // send their numbers to processors in the first row.
    if (y == 0 && x != 0)
    {
        MPI_Send(col_vector, nlocal_x, MPI_DOUBLE, reverse_rank, 111, comm);
    }


    if (x == 0 && y != 0)
    {
        MPI_Recv(row_vector, nlocal_y, MPI_DOUBLE, reverse_rank, 111, comm, &status);
    }

    if (x == 0 && y == 0)
    {
        for(int i= 0;i < nlocal_y; i++)
        {
            row_vector[i] = col_vector[i];
        }
    }


    // processors in the first row broadcast their elements to other processors in the same column
    int col_rank;

    MPI_Comm col_comm;
    MPI_Comm_split(comm, my2dcoords[1],my2drank,&col_comm);
    MPI_Comm_rank(col_comm, &col_rank);

    MPI_Bcast(row_vector, nlocal_y, MPI_DOUBLE, 0, col_comm);
    MPI_Comm_free(&col_comm);


}


void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm)
{
    // TODO

    int npes;
    int my2drank, my2dcoords[2];


    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &my2drank);
    MPI_Cart_coords(comm, my2drank, 2, my2dcoords);

    int x = my2dcoords[0];
    int y = my2dcoords[1];

    int q = sqrt(npes);

    int nlocal_x, nlocal_y;
    nlocal_x = block_decompose(n,q,x);
    nlocal_y = block_decompose(n,q,y);


    double local_x_distributed[nlocal_y];
    transpose_bcast_vector(n, local_x, local_x_distributed, comm);



    double local_y_distributed[nlocal_x];
    matrix_vector_mult(nlocal_x, nlocal_y, local_A, local_x_distributed, local_y_distributed);


    int rank;

    MPI_Comm row_comm;
    MPI_Comm_split(comm, my2dcoords[0],my2drank,&row_comm);
    MPI_Comm_rank(row_comm, &rank);

    int row_npes;
    MPI_Comm_size(row_comm, &row_npes);

    MPI_Reduce(local_y_distributed, local_y, nlocal_x, MPI_DOUBLE, MPI_SUM, 0, row_comm);
    MPI_Comm_free(&row_comm);



}

// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x,
                MPI_Comm comm, int max_iter, double l2_termination)
{
    // TODO
    
    
    int npes;
    int my2drank, my2dcoords[2];


    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &my2drank);
    MPI_Cart_coords(comm, my2drank, 2, my2dcoords);
    MPI_Status status;


    int x = my2dcoords[0];
    int y = my2dcoords[1];

    int q = sqrt(npes);

    int nlocal_x;
    nlocal_x = block_decompose(n,q,x);

    double* local_D = NULL;
    if (x == y || y == 0)
    {
        local_D = new double[nlocal_x];
    }
    

    // int nlocal = nlocal_x*nlocal_y; 

    if (x == y)
        {
            
            for (int i = 0; i < nlocal_x; i ++)
            {
                local_D[i] = local_A[i*nlocal_x + i];
            }

            if (x != 0)

            {
                int receive_rank;
                int receive_coord[2] = {x, 0};
    
                MPI_Cart_rank(comm, receive_coord, &receive_rank);
    
                MPI_Send(local_D, nlocal_x, MPI_DOUBLE, receive_rank, 111, comm);
            }
        }



    if (y == 0 && x != 0)
    {
       
        int send_rank;
        int send_coord[2] = {x, x};

        MPI_Cart_rank(comm, send_coord, &send_rank);
        MPI_Recv(local_D, nlocal_x, MPI_DOUBLE, send_rank, 111, comm, &status);
    }

    


    if (y == 0)
    {
        for (int i = 0; i < nlocal_x; i++)
            {
                local_x[i] = 0;
            }
    }


    int col_rank;

    MPI_Comm col_comm;
    MPI_Comm_split(comm, my2dcoords[1],my2drank,&col_comm);
    MPI_Comm_rank(col_comm, &col_rank);

 
    double l2 = 0;
    double local_l2 = 0;
    double* w = NULL;

    if (y == 0)
        {
            w = new double[nlocal_x];
            
        }
    
    
    for (int j = 0; j < max_iter; j++)
    {

        distributed_matrix_vector_mult(n, local_A, local_x, w, comm);

        if (y == 0)
        {
            for (int i = 0; i < nlocal_x; i++)
            {
                local_x[i] = (local_b[i] - w[i])/local_D[i] + local_x[i];
            }

            


            local_l2 = 0;

            for (int i = 0; i < nlocal_x; i++)
            {
                local_l2 = local_l2 + (local_b[i] - w[i])*(local_b[i] - w[i]);
            }
        }

        



        MPI_Reduce(&local_l2, &l2, 1, MPI_DOUBLE, MPI_SUM, 0, col_comm);
        MPI_Bcast(&l2, 1, MPI_DOUBLE, 0, comm);
        




        if (sqrt(l2) <= l2_termination)
        {
            MPI_Comm_free(&col_comm);
            return;
        }

    }

}


// wraps the distributed matrix vector multiplication
void mpi_matrix_vector_mult(const int n, double* A,
                            double* x, double* y, MPI_Comm comm)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_x = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &x[0], &local_x, comm);

    // allocate local result space
    double* local_y = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);

    // gather results back to rank 0
    gather_vector(n, local_y, y, comm);
}

// wraps the distributed jacobi function
void mpi_jacobi(const int n, double* A, double* b, double* x, MPI_Comm comm,
                int max_iter, double l2_termination)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_b = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &b[0], &local_b, comm);

    // allocate local result space
    double* local_x = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_jacobi(n, local_A, local_b, local_x, comm, max_iter, l2_termination);

    // gather results back to rank 0
    gather_vector(n, local_x, x, comm);
}
