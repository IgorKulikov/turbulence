#ifndef _save_hpp_
#define _save_hpp_
/*******************************************************************/
#include <stdio.h>
#include "global.hpp"
#include "system.hpp"
#include "parallel.hpp"

void save(real *a, real *bx, real *by, real timer)
{
	integer i, k, l;
	real total_density_a, x, y;
	FILE *foutxy;
	char filename[32];

	// create files
	sprintf(filename,"TXY_%3.6lf_p%d.dat",timer,rank);
	foutxy = fopen(filename,"w");
		
	// xy column density
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
		{
			x = (i+i_start_index)*h - h/2.0 - xm/2.0;
			y = k*h - h/2.0 - ym/2.0;

			total_density_a = a[index(i,k,NZ/2)];
			
			fprintf(foutxy,"%lf %lf %lf %lf %lf\n",
				x,y,total_density_a,bx[index(i,k,NZ/2)],by[index(i,k,NZ/2)]);
		}
			
	fclose(foutxy);
}

void save_log(FILE *fout, real timer)
{
	integer i, k, l;
	real x, y, z;
	real dmass, dlvx, dlvy, dlvz, drvx, drvy, drvz, denergy, dvelq, dkin, dint, dgrav;
	real gldmass, gldlvx, gldlvy, gldlvz, gldrvx, gldrvy, gldrvz, gldenergy;
	real gldkin, gldint, gldgrav;
	
	dmass = 0.0;
	dlvx  = 0.0;
	dlvy  = 0.0;
	dlvz  = 0.0;
	drvx  = 0.0;
	drvy  = 0.0;
	drvz  = 0.0;
	denergy = 0.0;
	dkin  = 0.0;
	dint  = 0.0;
	dgrav = 0.0;

	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				x = i*h - h/2.0 - xm/2.0;
				y = k*h - h/2.0 - ym/2.0;
				z = l*h - h/2.0 - zm/2.0;
				dvelq = Vx[index(i,k,l)] * Vx[index(i,k,l)] +
					    Vy[index(i,k,l)] * Vy[index(i,k,l)] +
					    Vz[index(i,k,l)] * Vz[index(i,k,l)];
				dmass   += R[index(i,k,l)];
				drvx    += RUx[index(i,k,l)];
				drvy    += RUy[index(i,k,l)];
				drvz    += RUz[index(i,k,l)];
				dlvx    += (y*RUz[index(i,k,l)] - z*RUy[index(i,k,l)]);
				dlvy    += (z*RUx[index(i,k,l)] - x*RUz[index(i,k,l)]);
				dlvz    += (x*RUy[index(i,k,l)] - y*RUx[index(i,k,l)]);
				dkin    += 0.5*R[index(i,k,l)]*dvelq;
				dint    += P[index(i,k,l)]/(GAMMA - 1.0);
				dgrav   += 0.5*R[index(i,k,l)]*Fi[index(i,k,l)];
				denergy += 0.5*R[index(i,k,l)]*dvelq + 
					P[index(i,k,l)]/(GAMMA - 1.0) + 
					0.5*R[index(i,k,l)]*Fi[index(i,k,l)] + 
					0.5 * (	Bx[index(i,k,l)] * Bx[index(i,k,l)] + 
							By[index(i,k,l)] * By[index(i,k,l)] +
							Bz[index(i,k,l)] * Bz[index(i,k,l)]);
			}
	dmass   *= h*h*h;
	drvx    *= h*h*h;
	drvy    *= h*h*h;
	drvz    *= h*h*h;
	dlvx    *= h*h*h;
	dlvy    *= h*h*h;
	dlvz    *= h*h*h;
	denergy *= h*h*h;
	dkin    *= h*h*h;
	dint    *= h*h*h;
	dgrav   *= h*h*h;

	MPI_Allreduce(&dmass,&gldmass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&drvx,&gldrvx,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&drvy,&gldrvy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&drvz,&gldrvz,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&dlvx,&gldlvx,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&dlvy,&gldlvy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&dlvz,&gldlvz,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&denergy,&gldenergy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&dkin,&gldkin,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&dint,&gldint,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&dgrav,&gldgrav,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	if(rank == 0)
	{
		fprintf(fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
				timer,gldmass,gldrvx,gldrvy,gldrvz,gldlvx,gldlvy,gldlvz,gldenergy,
				gldkin, gldint, gldgrav);
		fflush(fout);
	}

}

/*******************************************************************/
#endif