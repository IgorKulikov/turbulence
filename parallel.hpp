#ifndef _parallel_hpp_
#define _parallel_hpp_
/*******************************************************************/
#include "system.hpp"
#include "global.hpp"

// Tag of communiactions
#define TAGTOLEFT  1000
#define TAGTORIGHT 1001

// Number first cell
int i_start_index;

// Number of process and size of topology
int rank, size;		

// Send-Recv buffer
real *buffer;

/*********** Create forward plan ***********/
fftwnd_mpi_plan create_forward_plan(int iDim)
{
	return fftw3d_mpi_create_plan(MPI_COMM_WORLD, iDim, iDim, iDim,FFTW_FORWARD, FFTW_ESTIMATE);
}

/*********** Create backward plan ***********/
fftwnd_mpi_plan create_backward_plan(int iDim)
{
	return fftw3d_mpi_create_plan(MPI_COMM_WORLD, iDim, iDim, iDim,FFTW_BACKWARD, FFTW_ESTIMATE);
}

/*********** Destroy plan ***********/
void destroy_plan(fftwnd_mpi_plan plan)
{
	fftwnd_mpi_destroy_plan(plan);
}

// Parallel initial
void parallel_init()
{	
	// Total size for FFTW
	int iTotalSize;
	
	// Temporary variables for FFTW
	int local_ny_after_transpose,local_y_start_after_transpose;

	// Create plan
	plan  = create_forward_plan(FULL_NX);
	iplan = create_backward_plan(FULL_NX);
	fftwnd_mpi_local_sizes(plan, &locN, &i_start_index,
		&local_ny_after_transpose,&local_y_start_after_transpose,&iTotalSize);

	data = (fftw_complex*) malloc(sizeof(fftw_complex) * iTotalSize);
	schemas = (real*) malloc(sizeof(real) * FULL_NX);
	RFFTW = (real*) malloc(sizeof(real) * FULL_NX * FULL_NX * locN);
	FiFFTW = (real*) malloc(sizeof(real) * FULL_NX * FULL_NX * locN);

	NX = locN + 2;
	NY = FULL_NY + 2;
	NZ = FULL_NZ + 2;
	
	// Create buffer for communications
	buffer = new real[16*NY*NZ];
}

// Parallel destroy
void parallel_destroy()
{
	delete buffer;
}

/*******************************************************************/
#endif