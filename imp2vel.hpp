#ifndef _imp2vel_hpp_
#define _imp2vel_hpp_
/*******************************************************************/
#include "global.hpp"
#include "boundary.hpp"

real inv27p(real *r, integer i, integer k, integer l)
{
	return	8.0 * r[index(i,k,l)] +
			4.0 * ( r[index(i+1,k,l)] + r[index(i-1,k,l)] +
					r[index(i,k+1,l)] + r[index(i,k-1,l)] + 
					r[index(i,k,l+1)] + r[index(i,k,l-1)]		) +
			2.0 * ( r[index(i+1,k+1,l)] + r[index(i-1,k+1,l)] +
					r[index(i+1,k-1,l)] + r[index(i-1,k-1,l)] +
					r[index(i,k+1,l+1)] + r[index(i,k-1,l+1)] + 
					r[index(i,k+1,l-1)] + r[index(i,k-1,l-1)] + 
					r[index(i+1,k,l+1)] + r[index(i-1,k,l+1)] + 
					r[index(i+1,k,l-1)] + r[index(i-1,k,l-1)]	) +
			r[index(i+1,k+1,l+1)] + r[index(i-1,k+1,l+1)] + 
			r[index(i+1,k+1,l-1)] + r[index(i-1,k+1,l-1)] +
			r[index(i+1,k-1,l+1)] + r[index(i-1,k-1,l+1)] +
			r[index(i+1,k-1,l-1)] + r[index(i-1,k-1,l-1)] + EPSILON;
}


void imp2vel(real *r, real *rvx, real *rvy, real *rvz, real *vx, real *vy, real *vz)
{
	integer i, k, l;
	real density, dsummro, dinv;
	
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				density = r[index(i,k,l)];
				dsummro = inv27p(r,i,k,l);
				dinv = 4096.0 * density/dsummro/dsummro;
				vx[index(i,k,l)] = rvx[index(i,k,l)]*dinv;
				vy[index(i,k,l)] = rvy[index(i,k,l)]*dinv;	
				vz[index(i,k,l)] = rvz[index(i,k,l)]*dinv;	
			}

	boundary(vx);
	boundary(vy);
	boundary(vz);
}

void div2rho(real *ain, real *r, real *aout)
{
	integer i, k, l;
	real density, dsummro, dinv;
	
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				density = r[index(i,k,l)];
				dsummro = inv27p(r,i,k,l);
				dinv = 4096.0 * density/dsummro/dsummro;
				aout[index(i,k,l)] = ain[index(i,k,l)]*dinv;
			}

	boundary(aout);
}

/*******************************************************************/
#endif