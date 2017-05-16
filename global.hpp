#ifndef _global_hpp_
#define _global_hpp_
/*******************************************************************/
#include "system.hpp"
#include "mesh.hpp"

// Size of global mesh
#define FULL_NX 256
#define FULL_NY FULL_NX
#define FULL_NZ FULL_NX

// Courant–Friedrichs–Lewy condition
#define CFL 0.2

// Domain
real xm = 2.56, ym = 2.56, zm = 2.56;
real h = xm/FULL_NX;
real tau;

// Maximal time
real max_time = 1.5;

// Printer with "max_time = print_tau * iTimeCount"
real print_tau = 0.1;
integer iTimeCount = 15;

// Adiabatic index
#define GAMMA (5.0/3.0)

// Size of local mesh
integer NX, NY, NZ;

// Clear x-size of local mesh <<locN = NX -2>>
integer locN;

// Arrays for Hydrodynamics
real *R, *R_Next,						// density arrays
	 *RUx, *RUy, *RUz,					// impulses arrays
	 *RUx_Next, *RUy_Next, *RUz_Next,	 
	 *Vx, *Vy, *Vz,						// velocity arrays
	 *P, *P_Next,						// pressure arrays
	 *RE, *RE_Next,						// total energy arrays
	 *Bx, *By, *Bz,						// magnetic field
	 *Bx_Next, *By_Next, *Bz_Next,	
	 *Fi;								// gravity potential

// Arrays for mesh
mesh *parvx, *parvy, *parvz, *parp, *parbx, *parby, *parbz;

// Poisson arrays
real *RFFTW, *FiFFTW, *schemas;

// Plan for FFTW
fftwnd_mpi_plan plan, iplan;
fftw_complex *data;

// Index function
integer index(integer i, integer k, integer l)
{
	return i*NZ*NY+k*NZ+l;
}

/*******************************************************************/
#endif