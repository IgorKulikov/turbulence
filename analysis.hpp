#ifndef _analysis_hpp_
#define _analysis_hpp_
/*******************************************************************/
#include <stdio.h>
#include "global.hpp"
#include "system.hpp"
#include "parallel.hpp"

void analysis(real *r, real *vx, real *vy, real *vz, real *bx, real *by, real *bz, real timer)
{
	integer i, k, l;
	real mycosinus, malfven;
	FILE *foutxy;
	char filename[32];

	// create files
	sprintf(filename,"Tanalysis_%3.6lf_p%d.dat",timer,rank);
	foutxy = fopen(filename,"w");

	for(i=1;i<NX-1;i++)
    	    for(k=1;k<NY-1;k++)
		for(l=1;l<NZ-1;l++)
		{
		    mycosinus = ( vx[index(i,k,l)]*bx[index(i,k,l)] + 
				  vy[index(i,k,l)]*by[index(i,k,l)] + 
				  vz[index(i,k,l)]*bz[index(i,k,l)] ) / 
				  sqrt(	vx[index(i,k,l)]*vx[index(i,k,l)] + 
					vy[index(i,k,l)]*vy[index(i,k,l)] + 
					vz[index(i,k,l)]*vz[index(i,k,l)]) / 
				  sqrt(	bx[index(i,k,l)]*bx[index(i,k,l)] + 
					by[index(i,k,l)]*by[index(i,k,l)] + 
					bz[index(i,k,l)]*bz[index(i,k,l)]);
		    malfven = sqrt(4.0*PI*r[index(i,k,l)]) * 
				sqrt(	vx[index(i,k,l)]*vx[index(i,k,l)] + 
					vy[index(i,k,l)]*vy[index(i,k,l)] + 
					vz[index(i,k,l)]*vz[index(i,k,l)]) /
				sqrt(	bx[index(i,k,l)]*bx[index(i,k,l)] + 
					by[index(i,k,l)]*by[index(i,k,l)] + 
					bz[index(i,k,l)]*bz[index(i,k,l)]);
		
		    fprintf(foutxy,"%lf %lf %lf\n",
			r[index(i,k,l)],mycosinus,malfven);
		}

	fclose(foutxy);
}

/*******************************************************************/
#endif