#ifndef _problem_hpp_
#define _problem_hpp_
/*******************************************************************/
#include <math.h>
#include "parabola.hpp"
#include "system.hpp"
#include "global.hpp"
#include "parallel.hpp"
#include "boundary.hpp"

void load_problem()
{
	integer i, k, l;
	real x, y, z, rad, vteta, vphi, bteta, bphi, vel0, bfield0, dens0, press0;
	srand(73*rank+37);
	// init
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				x = (i+i_start_index)*h - h/2.0 - xm/2.0;
				y = k*h - h/2.0 - ym/2.0;
				z = l*h - h/2.0 - zm/2.0;
				
				vel0    = 1.0;
				bfield0 = 0.2;
				dens0   = 1.0;
				press0  = 1.0;

				vteta = random(0.0,2.0*PI);
				vphi  = random(0.0,2.0*PI);
				bteta = random(0.0,2.0*PI);
				bphi  = random(0.0,2.0*PI);

				R[index(i,k,l)]   = dens0;
				Vx[index(i,k,l)]  = vel0 * sin(vteta) * cos(vphi);
				Vy[index(i,k,l)]  = vel0 * sin(vteta) * sin(vphi);
				Vz[index(i,k,l)]  = vel0 * cos(vteta);
				Bx[index(i,k,l)]  = 0.0;
				By[index(i,k,l)]  = 0.0;
				Bz[index(i,k,l)]  = bfield0;
				P[index(i,k,l)]   = press0;

				RUx[index(i,k,l)] = R[index(i,k,l)] * Vx[index(i,k,l)];
				RUy[index(i,k,l)] = R[index(i,k,l)] * Vy[index(i,k,l)];
				RUz[index(i,k,l)] = R[index(i,k,l)] * Vz[index(i,k,l)];
				
				RE[index(i,k,l)] = P[index(i,k,l)]/(GAMMA - 1.0) + 
					0.5 * R[index(i,k,l)] * (	Vx[index(i,k,l)] * Vx[index(i,k,l)] + 
												Vy[index(i,k,l)] * Vy[index(i,k,l)] +
												Vz[index(i,k,l)] * Vz[index(i,k,l)] ) +
					0.5 * (	Bx[index(i,k,l)] * Bx[index(i,k,l)] + 
							By[index(i,k,l)] * By[index(i,k,l)] +
							Bz[index(i,k,l)] * Bz[index(i,k,l)]);
			}

	boundary(R);
	boundary(P);
	boundary(Vx);
	boundary(Vy);
	boundary(Vz);
	boundary(RUx);
	boundary(RUy);
	boundary(RUz);
	boundary(Bx);
	boundary(By);
	boundary(Bz);
	boundary(RE);	
}

/*******************************************************************/
#endif