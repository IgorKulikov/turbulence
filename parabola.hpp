#ifndef _parabola_hpp_
#define _parabola_hpp_
/*******************************************************************/
#include "system.hpp"

// type of parabola 
struct parabola
{	
	// parameters of parabola
	real ulr, ull, dq, dq6, ddq;

	// monotone
	void monotone(real uvalue, real uplus, real uminus, 
				  real dqplus, real dqminus, real h)
	{
		real eta, eta1 = 20.0, eta2 = 0.05, etas, ullp, ulrp;
		
		ullp = uvalue - 0.25 * dq;	
		ulrp = uvalue + 0.25 * dq;
		etas = -h * h * (dqplus*dqplus - dqminus*dqminus)/max(uplus-uminus,EPSILON);

		if( dqplus*dqminus > 0 ) etas = 0.0;
		if(	fabs(uplus-uminus) - 0.01*min(fabs(uplus),fabs(uminus),
										fabs(uplus)+fabs(uminus)) <= 0.0) etas = 0.0;

		eta = max(min(eta1*(etas-eta2),1.0,1.0),0.0);
		ull = (1.0 - eta) * ullp + eta * ull;
		ulr = (1.0 - eta) * ulrp + eta * ulr;
	}
	
	// right integral of parabola
	real right_parabola(real ksi, real h)
	{	
		return (ulr - 0.5 * ksi/h * (ddq - (1.0 - 2.0 * ksi / 3.0 / h) * dq6 ));	
	}

	// left integral of parabola
	real left_parabola(real ksi, real h)
	{	
		return (ull + 0.5 * ksi/h * (ddq + (1.0 - 2.0 * ksi / 3.0 / h) * dq6 ));	
	}

	// build local parabola
	void local_parabola(real uvalue)
	{
		ddq = ulr - ull;
		dq6 = 6.0 * (uvalue - (ulr + ull)/2.0);
	}

	// reconstruct parabola
	void reconstruct_parabola(real uvalue)
	{
		if( (ull - uvalue) * (uvalue - ulr) <= 0.0 )
		{
			ull = uvalue;
			ulr = uvalue;
		}
		else
		{
			if( ddq * dq6 > ddq * ddq )	
				ull = 3.0 * uvalue - 2.0 * ulr;
			else
				if( ddq * dq6 < -ddq * ddq )
					ulr = 3.0 * uvalue - 2.0 * ull;
		}
	}

};

/*******************************************************************/
#endif