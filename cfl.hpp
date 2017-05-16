#ifndef _cfl_hpp_
#define _cfl_hpp_
/*******************************************************************/
#include <math.h>
#include "global.hpp"
#include "system.hpp"

real computational_tau(real timer, real tend)
{
	integer i, k, l;
	real maxvcur, b2, c, ca, cf, cs, glmaxv;
	maxvcur = 0.0;

	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				b2 = (	Bx[index(i,k,l)]*Bx[index(i,k,l)] + 
						By[index(i,k,l)]*By[index(i,k,l)] +
						Bz[index(i,k,l)]*Bz[index(i,k,l)])/R[index(i,k,l)];
				c  = sqrt(GAMMA*P[index(i,k,l)]/R[index(i,k,l)]);
				ca = fabs(Bx[index(i,k,l)]/sqrt(R[i])) + 
					 fabs(By[index(i,k,l)]/sqrt(R[i])) + 
					 fabs(Bz[index(i,k,l)]/sqrt(R[i]));
				cf = sqrt(((c*c+b2) + sqrt((c*c+b2)*(c*c+b2) - 4*c*c*ca*ca))/2);
				cs = sqrt(((c*c+b2) - sqrt((c*c+b2)*(c*c+b2) - 4*c*c*ca*ca))/2);
				if( maxvcur < fabs(Vx[index(i,k,l)]) + c ) 
						maxvcur = fabs(Vx[index(i,k,l)]) + c;
				if( maxvcur < fabs(Vy[index(i,k,l)]) + c ) 
						maxvcur = fabs(Vy[index(i,k,l)]) + c;
				if( maxvcur < fabs(Vz[index(i,k,l)]) + c ) 
						maxvcur = fabs(Vz[index(i,k,l)]) + c;
				if( maxvcur < cf ) maxvcur = cf;
				if( maxvcur < cs ) maxvcur = cs;
			}
	MPI_Allreduce(&maxvcur,&glmaxv,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
	tau = CFL * h / glmaxv;
	if( timer + tau >= tend ) return tend - timer;
	else return tau;
}

/*******************************************************************/
#endif