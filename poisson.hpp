#ifndef _poisson_hpp_
#define _poisson_hpp_
/*******************************************************************/
#include "global.hpp"
#include "boundary.hpp"
#include "system.hpp"

/****************************************************************	
		Poisson solver 
****************************************************************/
void poisson_solver(real *R, real *Fi, fftw_complex *data, real *schemas,
			 fftwnd_mpi_plan plan, fftwnd_mpi_plan iplan, 
			 int iDim, real h, real x0)
{
	int local_nx, 
		local_x_start, 
		local_ny_after_transpose,
		local_y_start_after_transpose, 
		total_local_size;
	int x,y,z;
	int size, rank;
	real coeff;
	real dmove;
	real massa, global_massa;
	fftw_complex *cdata;
    
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);    

	/* get parameter distribution array */
	fftwnd_mpi_local_sizes(plan, &local_nx, &local_x_start,
                           &local_ny_after_transpose,
                           &local_y_start_after_transpose,
             			   &total_local_size);

        /* initialize schemas */
        for(x = 0 ; x < iDim ; x++ ) 
    	    schemas[x] = 1.0-(2.0/3.0) * (sin(x*PI/iDim)*sin(x*PI/iDim));

	/* set density */
	massa = 0.0;
	for (x = 0 ; x < local_nx ; x++)
	 for (y = 0 ; y < iDim ; y++)
	  for (z = 0 ; z < iDim ; z++)
		{
			data[x*iDim*iDim + y*iDim + z].re = R[x*iDim*iDim + y*iDim + z];
			data[x*iDim*iDim + y*iDim + z].im = 0.0;
			massa += R[x*iDim*iDim + y*iDim + z]*h*h*h;
		}		
	MPI_Allreduce(&massa,&global_massa,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	massa = global_massa/(12.56637061435917295385057);

	/* solve poisson equation */
	fftwnd_mpi(plan, 1, data, NULL, FFTW_TRANSPOSED_ORDER);
	cdata = (fftw_complex*) data;
	for (y = 0; y < local_ny_after_transpose; ++y)
    	    for (x = 0; x < iDim; ++x)
    		for (z = 0; z < iDim; ++z)
		{
			if( y+local_y_start_after_transpose+x+z == 0 ) 
				coeff = 0.0;
			else 
				coeff = -h*h/(6.0*(1.0 - schemas[x]*schemas[y+local_y_start_after_transpose]*schemas[z]));

			cdata[y*iDim*iDim + x*iDim + z].im *= coeff;
			cdata[y*iDim*iDim + x*iDim + z].re *= coeff;
		}
	fftwnd_mpi(iplan, 1, data, NULL, FFTW_TRANSPOSED_ORDER);


	for (x = 0; x < local_nx; x++)
	    for (y = 0; y < iDim; y++)
    		for (z = 0; z < iDim; z++)
			Fi[x*iDim*iDim + y*iDim + z] = 
				data[x*iDim*iDim + y*iDim + z].re/(iDim*iDim*iDim) + dmove;
			
}

// Laplass
void laplass()
{
	int i, k, l, indexl = 0;
	
	// Инициализация Пуассона
	indexl = 0;
	for(i=1;i<NX-1;i++)
	    for(k=1;k<NY-1;k++)
		for(l=1;l<NZ-1;l++)
		{
		    RFFTW[indexl] = 12.56637061435917295385057 * (R[index(i,k,l)] - 1.0);
		    indexl++;
		}

	// Счет уравнения Пуассона
	poisson_solver(RFFTW,FiFFTW,data,schemas,plan,iplan,FULL_NX,h,xm/2.0);
	
	// Сохранения уранвения Пуассона
	indexl = 0;
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				Fi[index(i,k,l)] = FiFFTW[indexl];
				indexl++;
			}
	boundary(Fi);
}

/*******************************************************************/
#endif