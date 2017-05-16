#ifndef _entropy_hpp_
#define _entropy_hpp_
/*******************************************************************/
#include "global.hpp"
#include "boundary.hpp"

void entropy(real *vx, real *vy, real *vz, 
	     real *bx, real *by, real *bz,
	     real *p,  real *re, real *r)
{
	integer i, k, l;
	real velq, mfld, trypres;
	
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				velq =	vx[index(i,k,l)]*vx[index(i,k,l)] + 
					vy[index(i,k,l)]*vy[index(i,k,l)] +
					vz[index(i,k,l)]*vz[index(i,k,l)];
				mfld = 0.5*(	bx[index(i,k,l)]*bx[index(i,k,l)] + 
						by[index(i,k,l)]*by[index(i,k,l)] +
						bz[index(i,k,l)]*bz[index(i,k,l)]);
				trypres = (GAMMA - 1.0) * 
					( re[index(i,k,l)] - 0.5 * r[index(i,k,l)] * velq - mfld );
				p[index(i,k,l)] = max(trypres,EPSILON);
			}

	boundary(p);
}

/*******************************************************************/
#endif