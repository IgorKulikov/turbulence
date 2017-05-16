#ifndef _system_hpp_
#define _system_hpp_
/*******************************************************************/
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <fftw.h>
#include <fftw_mpi.h>

// type of float 
#define real double

// type of integer
#define integer int

// PI value
#define PI 3.141592653589793238462

// Epsilon value
#define EPSILON 10e-5

// sqrt(4 \Pi) value
#define SQRT4PI			3.5449077018110320545963349666823

// 1/sqrt(2) value
#define ONEDIVSQRT2		0.70710678118654752440084436210485

// random function
real random(real tmin, real tmax)
{
	real h = tmax - tmin;
	return tmin + ((real)rand())/((real)RAND_MAX) * h;
}

// max function
real max(real t1, real t2)
{
	if(t1 > t2) return t1;
	return t2;
}

// min function
real min(real t1, real t2, real t3)
{
	real minimum = t1;
	if( t2 < minimum ) minimum = t2;
	if( t3 < minimum ) minimum = t3;
	return minimum;
}

// signum function
real signum(real t)
{
	if(t > 0.0) return 1;
	if(t == 0.0) return 0;
	return -1;
}

/*******************************************************************/
#endif