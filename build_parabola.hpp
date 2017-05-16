#ifndef _build_parabola_hpp_
#define _build_parabola_hpp_
/*******************************************************************/
#include "system.hpp"
#include "parallel.hpp"
#include "global.hpp"
#include "mesh.hpp"
#include "boundary.hpp"

// Build parabola
real build_parabola_getdq(real dplus, real dminus, real dzero)
{
	real dzerodq = 0.5 * (dplus - dminus);
	if( (dplus-dzero)*(dzero-dminus) > 0.0 )
		return signum(dzerodq) * 
			min(fabs(dzerodq),2.0*fabs(dplus-dzero),2.0*fabs(dzero-dminus));
	else
		return 0.0;
}

real build_parabola_getmedian(real uleft, real uright, real dqleft, real dqright)
{
	return 0.5 * (uright + uleft) - (dqright - dqleft)/6.0;
}

void build_parabola(real *u, mesh *a)
{
	integer i, k, l; 

	// init parabola 
	for(i=0;i<NX;i++)
		for(k=0;k<NY;k++)
			for(l=0;l<NZ;l++)
				a[index(i,k,l)].u = u[index(i,k,l)];

	// get delta
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				a[index(i,k,l)].parx.dq = 
					build_parabola_getdq(a[index(i+1,k,l)].u,a[index(i-1,k,l)].u,a[index(i,k,l)].u);
				a[index(i,k,l)].pary.dq = 
					build_parabola_getdq(a[index(i,k+1,l)].u,a[index(i,k-1,l)].u,a[index(i,k,l)].u);
				a[index(i,k,l)].parz.dq = 
					build_parabola_getdq(a[index(i,k,l+1)].u,a[index(i,k,l-1)].u,a[index(i,k,l)].u);
			}

	// boundary 
	boundary_mesh(a);

	// get median value
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				a[index(i,k,l)].parx.ull = 
					build_parabola_getmedian(a[index(i-1,k,l)].u,a[index(i,k,l)].u,
					 a[index(i-1,k,l)].parx.dq,a[index(i,k,l)].parx.dq);
				a[index(i,k,l)].parx.ulr = 
					build_parabola_getmedian(a[index(i,k,l)].u,a[index(i+1,k,l)].u,
					 a[index(i,k,l)].parx.dq,a[index(i+1,k,l)].parx.dq);

				a[index(i,k,l)].pary.ull = 
					build_parabola_getmedian(a[index(i,k-1,l)].u,a[index(i,k,l)].u,
					 a[index(i,k-1,l)].pary.dq,a[index(i,k,l)].pary.dq);
				a[index(i,k,l)].pary.ulr = 
					build_parabola_getmedian(a[index(i,k,l)].u,a[index(i,k+1,l)].u,
					 a[index(i,k,l)].pary.dq,a[index(i,k+1,l)].pary.dq);

				a[index(i,k,l)].parz.ull = 
					build_parabola_getmedian(a[index(i,k,l-1)].u,a[index(i,k,l)].u,
					 a[index(i,k,l-1)].parz.dq,a[index(i,k,l)].parz.dq);
				a[index(i,k,l)].parz.ulr = 
					build_parabola_getmedian(a[index(i,k,l)].u,a[index(i,k,l+1)].u,
					 a[index(i,k,l)].parz.dq,a[index(i,k,l+1)].parz.dq);
			}

	// get first parabola
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				a[index(i,k,l)].parx.local_parabola(a[index(i,k,l)].u);
				a[index(i,k,l)].pary.local_parabola(a[index(i,k,l)].u);
				a[index(i,k,l)].parz.local_parabola(a[index(i,k,l)].u);
			}

	// reconstruct parabola
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				a[index(i,k,l)].parx.reconstruct_parabola(a[index(i,k,l)].u);
				a[index(i,k,l)].pary.reconstruct_parabola(a[index(i,k,l)].u);
				a[index(i,k,l)].parz.reconstruct_parabola(a[index(i,k,l)].u);
			}

	// monotone parabola
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				a[index(i,k,l)].parx.monotone(
					a[index(i,k,l)].u,a[index(i+1,k,l)].u,a[index(i-1,k,l)].u,
					a[index(i+1,k,l)].parx.dq,a[index(i-1,k,l)].parx.dq,h);

				a[index(i,k,l)].pary.monotone(
					a[index(i,k,l)].u,a[index(i,k+1,l)].u,a[index(i,k-1,l)].u,
					a[index(i,k+1,l)].pary.dq,a[index(i,k-1,l)].pary.dq,h);
				
				a[index(i,k,l)].parz.monotone(
					a[index(i,k,l)].u,a[index(i,k,l+1)].u,a[index(i,k,l-1)].u,
					a[index(i,k,l+1)].parz.dq,a[index(i,k,l-1)].parz.dq,h);				
			}

	// get second parabola
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				a[index(i,k,l)].parx.local_parabola(a[index(i,k,l)].u);
				a[index(i,k,l)].pary.local_parabola(a[index(i,k,l)].u);
				a[index(i,k,l)].parz.local_parabola(a[index(i,k,l)].u);
			}

	// boundary 
	boundary_mesh(a);
}

/*******************************************************************/
#endif