CSE 6220 High Performance Computing Final Project
=================================
## Project Overview
* Inplemented Jacobi's methold for iteratively solving large system of linear equations sequentially and paralelly using MPI in C++
* Tested the performance of the parallel implementation against the sequential implementation for various input size (max n = 1000), number of processors (max p = 64), and difficulties (max d = 0.9)

* Summarized and analyzied the results in the [report](https://github.com/sliao7/High-Performance-Computing/blob/master/jacobi_iterative/assign_3/report.pdf)

## Code Structure

All the code is located at the root level of the project.
The `gtest` folder contains the Google Test Unit Testing framework version 1.7.

There are multiple header and .cpp files, your implementation will go
into the following files:

- `jacobi.cpp`: Implement the sequential algorithm for Jacobi's method according
  to the function declarations in `jacobi.h`
- `mpi_jacobi.cpp`: Implement the parallel algorithm according to the
  declarations in `mpi_jacobi.h`.
- `utils.h` and `utils.cpp`: Implement common utility functions in these files
- `seq_tests.cpp`: Unit tests for the sequential code, implement your own
  test cases in here
- `mpi_tests.cpp`: Unit tests for the parallel MPI code. Implement your own
  test cases for all declared functions in here.


Other files containing code that you should not change are:

- `main.cpp`: Implements code for the main executable `jacobi`. This does
  input/output reading and calling of the actual functions.
- `io.h`: implements IO functions and random input generation
- `mpi_gtest.cpp`: MPI wrapper for the GTest framework


Utility scripts (you may play around with these to generate your own custom
input):

- `generate_input.py`: Python script to generate inputs. You can modify this
  code to generate different inputs.
- `check_output.py`: Checks whether the output from `jacobi` is correct, by
  comparing the output with pythons numpy implementation.


## Compiling

In order to compile everything, simply run
```sh
make all
```


## Running Tests

For running all tests do:
```sh
make test
```

You can also run the two tests (sequential/MPI) separately by either:
```sh
./seq_tests
```
or
```sh
mpirun -np 4 ./mpi_tests
```
