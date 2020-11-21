/**
 * @file    mpi_tests.cpp
 * @ingroup group
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   GTest Unit Tests for the parallel MPI code.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
/*
 * Add your own test cases here. We will test your final submission using
 * a more extensive tests suite. Make sure your code works for many different
 * input cases.
 *
 * Note:
 * The google test framework is configured, such that
 * only errors from the processor with rank = 0 are shown.
 */

#include <mpi.h>
#include <gtest/gtest.h>

#include <math.h>
#include "jacobi.h"
#include "mpi_jacobi.h"
#include "utils.h"
#include "io.h"

/**
 * @brief Creates and returns the square 2d grid communicator for MPI_COMM_WORLD
 */
void get_grid_comm(MPI_Comm* grid_comm)
{
    // get comm size and rank
    int rank, p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int q = (int)sqrt(p);
    ASSERT_EQ(q*q, p) << "Number of processors must be a perfect square.";

    // split into grid communicator
    int dims[2] = {q, q};
    int periods[2] = {0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, grid_comm);
}

// test parallel MPI matrix vector multiplication
TEST(MpiTest, distribute_vector)
{

//     // distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
//     // simple 4 by 4 input matrix
    double input_vector[5] =  {6., 25., -11., 15., 21};
    int n = 5;

    double* local_vector = NULL;
   
    

    // get grid communicator
    MPI_Comm grid_comm;
    get_grid_comm(&grid_comm);

    // int rank,p;
    // MPI_Comm_size(grid_comm,&p);
    // MPI_Comm_rank(grid_comm, &rank);

    distribute_vector(n, input_vector, &local_vector, grid_comm);    
    
}

TEST(MpiTest, gather_vector)
{

//     // distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
//     // simple 4 by 4 input matrix
    double input_vector[5] =  {6., 25., -11., 15., 21};
    int n = 5;

    double* local_vector = NULL;
    double output_vector[n];
    

    // get grid communicator
    MPI_Comm grid_comm;
    get_grid_comm(&grid_comm);

    // int rank,p;
    // MPI_Comm_size(grid_comm,&p);
    // MPI_Comm_rank(grid_comm, &rank);

    distribute_vector(n, input_vector, &local_vector, grid_comm);
    gather_vector(n, local_vector, output_vector, grid_comm);
    
    
}

TEST(MpiTest, distribute_matrix)
{

//     // distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
//     // simple 4 by 4 input matrix
    double input_matrix[25] =  {6., 25., -11., 15., 21.,
                                4., 3., 6., 5., 9., 
                               7., 12., 43., 23., 65.,
                               93., 29., 45., 22., 14.,
                               39., 59., 44., 30., 19.};
    int n = 5;

    double* local_matrix = NULL;
   
    

    // get grid communicator
    MPI_Comm grid_comm;
    get_grid_comm(&grid_comm);

 
    distribute_matrix(n, input_matrix, &local_matrix, grid_comm);    
    
}

TEST(MpiTest, transpose_bcast_vector)
{

//     // distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
//     // simple 4 by 4 input matrix
    // double vector[25] =  {6., 25., -11., 15., 21.,
    //                             4., 3., 6., 5., 9., 
    //                            7., 12., 43., 23., 65.,
    //                            93., 29., 45., 22., 14.,
    //                            39., 59., 44., 30., 19.};
    double vector[5] = {1., 2., 3., 4., 5.};
    int n = 5;
   
    

    // get grid communicator
    MPI_Comm grid_comm;
    get_grid_comm(&grid_comm);

    int my2drank, my2dcoords[2];
    int npes;

    MPI_Comm_size(grid_comm, &npes);
    MPI_Comm_rank(grid_comm, &my2drank);
    MPI_Cart_coords(grid_comm, my2drank, 2, my2dcoords);

    int x = my2dcoords[0];
    int y = my2dcoords[1];

    
    int nlocal_y;

    int q = sqrt(npes);
   
    nlocal_y = block_decompose(n,q,y);

    double* col_vector;

    int counts[q];


    for (int i = 0; i < q; i++)
    {       
        counts[i] = block_decompose(n, q, i);      
    }



    int displ[q];
    displ[0] = 0;
    for (int i = 1; i < q; i++)
    {
        displ[i] = displ[i-1] + counts[i-1];
        
    }


    if (y == 0)
    {
        col_vector = &vector[displ[x]];
    }
    else
    {
        col_vector = NULL;
    }
    


    double row_vector[nlocal_y];

    
    transpose_bcast_vector(n, col_vector, row_vector, grid_comm);
    
     
    
}

TEST(MpiTest, distributed_matrix_vector_mult)
{

//     // distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
//     // simple 4 by 4 input matrix
    double input_matrix[25] =  {6., 25., -11., 15., 21.,
                                4., 3., 6., 5., 9., 
                               7., 12., 43., 23., 65.,
                               93., 29., 45., 22., 14.,
                               39., 59., 44., 30., 19.};
    double input_x[5] = {1.0, 2.0, 3.0, 4., 5.};
    int n = 5;

    double* local_A = NULL;
    double* local_x = NULL;
   
    

    // get grid communicator
    MPI_Comm grid_comm;
    get_grid_comm(&grid_comm);


    distribute_matrix(n, input_matrix, &local_A, grid_comm);   
    distribute_vector(n, input_x, &local_x, grid_comm);

    int my2drank, my2dcoords[2];
    int npes;

    MPI_Comm_size(grid_comm, &npes);
    MPI_Comm_rank(grid_comm, &my2drank);
    MPI_Cart_coords(grid_comm, my2drank, 2, my2dcoords);

    int x = my2dcoords[0];
   

    
    int nlocal_x;

    int q = sqrt(npes);
    nlocal_x = block_decompose(n,q,x);
   

    double local_y[nlocal_x]; 
    distributed_matrix_vector_mult(n, local_A, local_x, local_y, grid_comm);



    
}






// test parallel MPI matrix vector multiplication
TEST(MpiTest, MatrixVectorMult1)
{
    // simple 4 by 4 input matrix
    double A[4*4] = {10., -1., 2., 0.,
                           -1., 11., -1., 3.,
                           2., -1., 10., -1.,
                           0.0, 3., -1., 8.};
    double x[4] =  {6., 25., -11., 15.};
    double y[4];
    double expected_y[4] = {13.,  325., -138.,  206.};
    int n = 4;

    // get grid communicator
    MPI_Comm grid_comm;
    get_grid_comm(&grid_comm);

    // testing sequential matrix multiplication
    mpi_matrix_vector_mult(n, A, x, y, grid_comm);

    // checking if all values are correct (up to some error value)
    for (int i = 0; i < n; ++i)
    {
        EXPECT_NEAR(expected_y[i], y[i], 1e-10) << " element y[" << i << "] is wrong";
    }
}


// // test parallel MPI matrix vector multiplication
TEST(MpiTest, Jacobi1)
{
    // simple 4 by 4 input matrix
    double A[4*4] = {10., -1., 2., 0.,
                           -1., 11., -1., 3.,
                           2., -1., 10., -1.,
                           0.0, 3., -1., 8.};
    double b[4] =  {6., 25., -11., 15.};
    double x[4];
    double expected_x[4] = {1.0,  2.0, -1.0, 1.0};
    int n = 4;

    // get grid communicator
    MPI_Comm grid_comm;
    get_grid_comm(&grid_comm);

    // testing sequential matrix multiplication
    mpi_jacobi(n, A, b, x, grid_comm);
   

    // checking if all values are correct (up to some error value)
    for (int i = 0; i < n; ++i)
    {
        EXPECT_NEAR(expected_x[i], x[i], 1e-5) << " element y[" << i << "] is wrong";
    }
}


// /**
//  * Test the parallel code and compare the results with the sequential code.
//  */
TEST(MpiTest, JacobiCrossTest1)
{
    // test random matrixes, test parallel code with sequential solutions
    std::vector<double> A;
    std::vector<double> b;
    std::vector<double> mpi_x;

    // get grid communicator
    MPI_Comm grid_comm;
    get_grid_comm(&grid_comm);
    int rank;
    MPI_Comm_rank(grid_comm, &rank);

    int n = 36;
    // initialize data only on rank 0
    if (rank == 0)
    {
        A = diag_dom_rand(n);
        b = randn(n, 100.0, 50.0);
    }

    // getting sequential results
    std::vector<double> x;
    if (rank == 0)
    {
        x.resize(n);
        jacobi(n, &A[0], &b[0], &x[0]);
    }

    // parallel jacobi
    if (rank == 0)
        mpi_x.resize(n);
    mpi_jacobi(n, &A[0], &b[0], &mpi_x[0], grid_comm);

    if (rank == 0)
    {
        // checking if all values are correct (up to some error value)
        for (int i = 0; i < n; ++i)
        {
            EXPECT_NEAR(x[i], mpi_x[i], 1e-8) << " MPI solution x[" << i << "] differs from sequential result";
        }
    }
}
