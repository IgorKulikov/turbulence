#ifndef _riemann_hpp_
#define _riemann_hpp_
/*******************************************************************/
#include <math.h>
#include "system.hpp"
#include "global.hpp"
#include "build_parabola.hpp"
#include "parabola.hpp"

void mhd_riemann(real rho_plus, real rho_minus, real p_plus, real p_minus,
				 real vx_plus, real vx_minus, real vy_plus, real vy_minus,
				 real vz_plus, real vz_minus, real bx_plus, real bx_minus,
				 real by_plus, real by_minus, real bz_plus, real bz_minus,
				 parabola pp_plus, parabola pp_minus, 
				 parabola pvx_plus, parabola pvx_minus,
				 parabola pvy_plus, parabola pvy_minus,
				 parabola pvz_plus, parabola pvz_minus,
				 parabola pby_plus, parabola pby_minus,
				 parabola pbz_plus, parabola pbz_minus,
				 real& fvx, real& fvy, real& fvz, real& fp, 
				 real& fbx, real& fby, real& fbz)
{
	real sqrp, sqrm, mulsqr, r, p, vx, vy, vz, bx, by, bz; 
	real cs, ca, cf, ct, b2, sgnBx, betay, betaz, alphaf, alphat;
	
	real VX_MCFT,VY_MCFT,VZ_MCFT,BY_MCFT,BZ_MCFT,P_MCFT,
		 VX_PCFT,VY_PCFT,VZ_PCFT,BY_PCFT,BZ_PCFT,P_PCFT,
		 VX_MCST,VY_MCST,VZ_MCST,BY_MCST,BZ_MCST,P_MCST,
		 VX_PCST,VY_PCST,VZ_PCST,BY_PCST,BZ_PCST,P_PCST,
		 VX_MCAT,VY_MCAT,VZ_MCAT,BY_MCAT,BZ_MCAT,P_MCAT,
		 VX_PCAT,VY_PCAT,VZ_PCAT,BY_PCAT,BZ_PCAT,P_PCAT;

	sqrp = sqrt(max(rho_plus,EPSILON));
	sqrm = sqrt(max(rho_minus,EPSILON));
	mulsqr = sqrp * sqrm;
	bx  = 0.5 * (bx_plus + bx_minus);
	sgnBx = signum(bx);

	r   = max((sqrp * rho_plus + sqrm * rho_minus)/(sqrp + sqrm),EPSILON);
	p   = max((sqrp * p_plus + sqrm * p_minus)/(sqrp + sqrm),EPSILON);
	vx  = (sqrp * vx_plus + sqrm * vx_minus)/(sqrp + sqrm);
	vy  = (sqrp * vy_plus + sqrm * vy_minus)/(sqrp + sqrm);
	vz  = (sqrp * vz_plus + sqrm * vz_minus)/(sqrp + sqrm);
	by  = (sqrm * by_plus + sqrp * by_minus)/(sqrp + sqrm);
	bz  = (sqrm * bz_plus + sqrp * bz_minus)/(sqrp + sqrm);

	b2 = (bx*bx + by*by + bz*bz)/r;
	cs = sqrt(GAMMA * p / r);
	ca = fabs(bx / sqrt(r));
	cf = sqrt(((cs*cs+b2) + sqrt((cs*cs+b2)*(cs*cs+b2) - 4.0*cs*cs*ca*ca))/2.0);
	ct = sqrt(((cs*cs+b2) - sqrt((cs*cs+b2)*(cs*cs+b2) - 4.0*cs*cs*ca*ca))/2.0);

	if( by*by + bz*bz > EPSILON )
	{
		betay = by/sqrt(by*by + bz*bz);
		betaz = bz/sqrt(by*by + bz*bz);
	}
	else
	{
		betay = ONEDIVSQRT2;
		betaz = ONEDIVSQRT2;
	}

	if( (by*by + bz*bz > EPSILON) || fabs(GAMMA*p - bx*bx) > EPSILON )
	{
		alphaf = sqrt(cs*cs-ct*ct + EPSILON)/sqrt(cf*cf-ct*ct);
		alphat = sqrt(cf*cf-cs*cs + EPSILON)/sqrt(cf*cf-ct*ct);
	}
	else
	{
		alphaf = ONEDIVSQRT2;
		alphat = ONEDIVSQRT2;
	}

	/* Analytical solution */
	VX_MCFT = pvx_minus.right_parabola(cf * tau, h);  
	VY_MCFT = pvy_minus.right_parabola(cf * tau, h);
	VZ_MCFT = pvz_minus.right_parabola(cf * tau, h);
	BY_MCFT = pby_minus.right_parabola(cf * tau, h);
	BZ_MCFT = pbz_minus.right_parabola(cf * tau, h);
	P_MCFT  =  pp_minus.right_parabola(cf * tau, h);

	VX_PCFT = pvx_plus.left_parabola(cf * tau, h);  
	VY_PCFT = pvy_plus.left_parabola(cf * tau, h);  
	VZ_PCFT = pvz_plus.left_parabola(cf * tau, h);  
	BY_PCFT = pby_plus.left_parabola(cf * tau, h);  
	BZ_PCFT = pbz_plus.left_parabola(cf * tau, h);  
	P_PCFT  =  pp_plus.left_parabola(cf * tau, h);  

	VX_MCST = pvx_minus.right_parabola(ct * tau, h);
	VY_MCST = pvy_minus.right_parabola(ct * tau, h);
	VZ_MCST = pvz_minus.right_parabola(ct * tau, h);
	BY_MCST = pby_minus.right_parabola(ct * tau, h);
	BZ_MCST = pbz_minus.right_parabola(ct * tau, h);
	P_MCST  =  pp_minus.right_parabola(ct * tau, h);

	VX_PCST = pvx_plus.left_parabola(ct * tau, h);  
	VY_PCST = pvy_plus.left_parabola(ct * tau, h);  
	VZ_PCST = pvz_plus.left_parabola(ct * tau, h);  
	BY_PCST = pby_plus.left_parabola(ct * tau, h);  
	BZ_PCST = pbz_plus.left_parabola(ct * tau, h);  
	P_PCST  =  pp_plus.left_parabola(ct * tau, h); 

	VX_MCAT = pvx_minus.right_parabola(ca * tau, h);
	VY_MCAT = pvy_minus.right_parabola(ca * tau, h);
	VZ_MCAT = pvz_minus.right_parabola(ca * tau, h);
	BY_MCAT = pby_minus.right_parabola(ca * tau, h);
	BZ_MCAT = pbz_minus.right_parabola(ca * tau, h);
	P_MCAT  =  pp_minus.right_parabola(ca * tau, h);

	VX_PCAT = pvx_plus.left_parabola(ca * tau, h);  
	VY_PCAT = pvy_plus.left_parabola(ca * tau, h);  
	VZ_PCAT = pvz_plus.left_parabola(ca * tau, h);  
	BY_PCAT = pby_plus.left_parabola(ca * tau, h);  
	BZ_PCAT = pbz_plus.left_parabola(ca * tau, h);  
	P_PCAT  =  pp_plus.left_parabola(ca * tau, h);  

	fvx = alphat * alphat * ct * ct / (cs *cs) * VX_MCST / 2.0 
		+ alphaf * alphaf * cf * cf / (cs *cs) * VX_PCFT / 2.0 
		- alphaf * cf * alphat * ct * betay * sgnBx / (cs *cs) * VY_MCFT / 2.0 
		- alphaf * cf * alphat * ct * betaz * sgnBx / (cs *cs) * VZ_MCFT / 2.0 
		+ alphaf * cf * alphat * betay / cs / sqrt(r) * BY_MCFT / 2.0 
		+ alphaf * cf * alphat * betaz / cs / sqrt(r) * BZ_MCFT / 2.0 
		+ alphaf * alphaf * cf / r / (cs *cs) * P_MCFT / 2.0 
		+ alphat * alphat * ct * ct / (cs *cs) * VX_PCST / 2.0 
		+ alphaf * cf * alphat * ct * betay * sgnBx / (cs *cs) * VY_PCST / 2.0 
		+ alphaf * cf * alphat * ct * betaz * sgnBx / (cs *cs) * VZ_PCST / 2.0 
		+ alphat * ct * alphaf * betay / cs / sqrt(r) * BY_PCST / 2.0 
		+ alphat * ct * alphaf * betaz / cs / sqrt(r) * BZ_PCST / 2.0 
		- alphat * alphat * ct / r / (cs *cs) * P_PCST / 2.0 
		+ alphaf * cf * alphat * ct * betaz * sgnBx / (cs *cs) * VZ_MCST / 2.0 
		- alphat * ct * alphaf * betay / cs / sqrt(r) * BY_MCST / 2.0 
		- alphat * ct * alphaf * betaz / cs / sqrt(r) * BZ_MCST / 2.0 
		+ alphat * alphat * ct / r / (cs *cs) * P_MCST / 2.0 
		+ alphaf * alphaf * cf * cf / (cs *cs) * VX_MCFT / 2.0 
		+ alphaf * cf * alphat * ct * betay * sgnBx / (cs *cs) * VY_MCST / 2.0 
		- alphaf * cf * alphat * ct * betay * sgnBx / (cs *cs) * VY_PCFT / 2.0 
		- alphaf * cf * alphat * ct * betaz * sgnBx / (cs *cs) * VZ_PCFT / 2.0 
		- alphaf * cf * alphat * betay / cs / sqrt(r) * BY_PCFT / 2.0 
		- alphaf * cf * alphat * betaz / cs / sqrt(r) * BZ_PCFT / 2.0 
		- alphaf * alphaf * cf * P_PCFT / r / (cs *cs) / 2.0;

	fvy = betaz * betaz * sgnBx * sgnBx * VY_MCAT / 2.0 
		- sgnBx * alphat * alphat * ct * betay * betaz / cs / sqrt(r) * BZ_MCFT / 2.0 
		+ betaz * betaz * sgnBx * sgnBx * VY_PCAT / 2.0 
		+ sgnBx * alphaf * cf * betay * alphat * ct / (cs *cs) * VX_MCST / 2.0 
		+ sgnBx * alphaf * alphaf * cf * betay * betay / cs / sqrt(r) * BY_PCST / 2.0 
		+ alphat * alphat * ct * ct * betay * betaz * sgnBx * sgnBx / (cs *cs) * VZ_MCFT / 2.0 
		+ sgnBx * alphaf * cf * betay * alphat / r / (cs *cs) * P_MCST / 2.0 
		- sgnBx * alphat * ct * betay * alphaf / r / (cs *cs) * P_MCFT / 2.0 
		+ alphaf * alphaf * cf * cf * betay * betaz * sgnBx * sgnBx / (cs *cs) * VZ_MCST / 2.0 
		+ sgnBx * betaz / sqrt(r) * betay * BZ_MCAT / 2.0 
		- sgnBx * alphaf * cf * betay * alphat * ct / (cs *cs) * VX_PCFT / 2.0 
		+ sgnBx * alphaf * alphaf * cf * betay * betaz / cs / sqrt(r) * BZ_PCST / 2.0 
		- sgnBx * sgnBx * betaz * betay * VZ_MCAT / 2.0 
		- sgnBx * alphat * alphat * ct * betay * betay / cs / sqrt(r) * BY_MCFT / 2.0 
		- sgnBx * betaz * betaz / sqrt(r) * BY_MCAT / 2.0 
		- sgnBx * alphaf * cf * betay * alphat * ct / (cs *cs) * VX_MCFT / 2.0 
		+ alphaf * alphaf * cf * cf * betay * betay * sgnBx * sgnBx / (cs *cs) * VY_MCST / 2.0 
		+ alphaf * alphaf * cf * cf * betay * betaz * sgnBx * sgnBx / (cs *cs) * VZ_PCST / 2.0 
		+ sgnBx * alphat * alphat * ct * betay * betay / cs / sqrt(r) * BY_PCFT / 2.0 
		+ sgnBx * alphaf * cf * betay * alphat * ct / (cs *cs) * VX_PCST / 2.0 
		+ alphat * alphat * ct * ct * betay * betay * sgnBx * sgnBx / (cs *cs) * VY_MCFT / 2.0 
		+ sgnBx * betaz * betaz / sqrt(r) * BY_PCAT / 2.0 
		- sgnBx * alphaf * alphaf * cf * betay * betaz / cs / sqrt(r) * BZ_MCST / 2.0 
		+ alphaf * alphaf * cf * cf * betay * betay * sgnBx * sgnBx / (cs *cs) * VY_PCST / 2.0 
		- sgnBx * alphaf * cf * betay * alphat / r / (cs *cs) * P_PCST / 2.0 
		- sgnBx * sgnBx * betaz * betay * VZ_PCAT / 2.0 
		+ alphat * alphat * ct * ct * betay * betaz * sgnBx * sgnBx / (cs *cs) * VZ_PCFT / 2.0 
		- sgnBx * alphaf * alphaf * cf * betay * betay / cs / sqrt(r) * BY_MCST / 2.0 
		- sgnBx * betaz / sqrt(r) * betay * BZ_PCAT / 2.0 
		+ sgnBx * alphat * ct * betay * alphaf * P_PCFT / r / (cs *cs) / 2.0 
		+ sgnBx * alphat * alphat * ct * betay * betaz / cs / sqrt(r) * BZ_PCFT / 2.0 
		+ alphat * alphat * ct * ct * betay * betay * sgnBx * sgnBx / (cs *cs) * VY_PCFT / 2.0;

	fvz = alphat * alphat * ct * ct * betaz * betaz * sgnBx * sgnBx / (cs *cs) * VZ_PCFT / 2.0 
		+ sgnBx * alphat * ct * betaz * alphaf * P_PCFT / r / (cs *cs) / 2.0 
		+ alphaf * alphaf * cf * cf * betaz * betaz * sgnBx * sgnBx / (cs *cs) * VZ_PCST / 2.0 
		+ sgnBx * alphaf * cf * betaz * alphat / r / (cs *cs) * P_MCST / 2.0 
		- sgnBx * alphaf * cf * betaz * alphat * ct / (cs *cs) * VX_PCFT / 2.0 
		+ sgnBx * alphat * alphat * ct * betay * betaz / cs / sqrt(r) * BY_PCFT / 2.0 
		+ sgnBx * alphaf * alphaf * cf * betay * betaz / cs / sqrt(r) * BY_PCST / 2.0 
		+ alphaf * alphaf * cf * cf * betay * betaz * sgnBx * sgnBx / (cs *cs) * VY_PCST / 2.0 
		- sgnBx * alphaf * cf * betaz * alphat / r / (cs *cs) * P_PCST / 2.0 
		+ sgnBx * betay * betay / sqrt(r) * BZ_PCAT / 2.0 
		- sgnBx * betay * betay / sqrt(r) * BZ_MCAT / 2.0 
		+ alphat * alphat * ct * ct * betay * betaz * sgnBx * sgnBx / (cs *cs) * VY_MCFT / 2.0 
		+ sgnBx * alphat * alphat * ct * betaz * betaz / cs / sqrt(r) * BZ_PCFT / 2.0 
		+ betay * betay * sgnBx * sgnBx * VZ_MCAT / 2.0 
		- sgnBx * alphaf * alphaf * cf * betaz * betaz / cs / sqrt(r) * BZ_MCST / 2.0 
		+ betay * betay * sgnBx * sgnBx * VZ_PCAT / 2.0 
		+ alphat * alphat * ct * ct * betay * betaz * sgnBx * sgnBx / (cs *cs) * VY_PCFT / 2.0 
		- sgnBx * sgnBx * betaz * betay * VY_PCAT / 2.0 
		- sgnBx * betaz / sqrt(r) * betay * BY_PCAT / 2.0 
		- sgnBx * alphaf * cf * betaz * alphat * ct / (cs *cs) * VX_MCFT / 2.0 
		+ alphat * alphat * ct * ct * betaz * betaz * sgnBx * sgnBx / (cs *cs) * VZ_MCFT / 2.0 
		+ alphaf * alphaf * cf * cf * betay * betaz * sgnBx * sgnBx / (cs *cs) * VY_MCST / 2.0 
		- sgnBx * sgnBx * betaz * betay * VY_MCAT / 2.0 
		- sgnBx * alphat * alphat * ct * betaz * betaz / cs / sqrt(r) * BZ_MCFT / 2.0 
		+ sgnBx * betaz / sqrt(r) * betay * BY_MCAT / 2.0 
		+ sgnBx * alphaf * alphaf * cf * betaz * betaz / cs / sqrt(r) * BZ_PCST / 2.0 
		- sgnBx * alphat * ct * betaz * alphaf / r / (cs *cs) * P_MCFT / 2.0 
		+ alphaf * alphaf * cf * cf * betaz * betaz * sgnBx * sgnBx / (cs *cs) * VZ_MCST / 2.0 
		+ sgnBx * alphaf * cf * betaz * alphat * ct / (cs *cs) * VX_MCST / 2.0 
		- sgnBx * alphat * alphat * ct * betay * betaz / cs / sqrt(r) * BY_MCFT / 2.0 
		- sgnBx * alphaf * alphaf * cf * betay * betaz / cs / sqrt(r) * BY_MCST / 2.0 
		+ sgnBx * alphaf * cf * betaz * alphat * ct / (cs *cs) * VX_PCST / 2.0;

	fby = alphaf * alphaf * betay * cf * betaz * sgnBx * sqrt(r) / cs * VZ_PCST / 2.0 
		+ alphat * alphat * betay * betaz * BZ_MCFT / 2.0 
		+ alphaf * alphaf * betay * betaz * BZ_PCST / 2.0 
		- betaz * betaz * sqrt(r) * sgnBx * VY_MCAT / 2.0 
		+ alphat * alphat * betay * betay * BY_MCFT / 2.0 
		- alphaf * betay * alphat * ct * sqrt(r) / cs * VX_MCST / 2.0 
		- alphat * betay * alphaf * cf * sqrt(r) / cs * VX_PCFT / 2.0 
		+ betaz * betaz * BY_MCAT / 2.0 
		+ betaz * betaz * BY_PCAT / 2.0 
		- alphat * betay * alphaf / sqrt(r) / cs * P_PCST / 2.0 
		+ betaz * sqrt(r) * betay * sgnBx * VZ_MCAT / 2.0 
		+ alphaf * alphaf * betay * betaz * BZ_MCST / 2.0 
		+ betaz * betaz * sqrt(r) * sgnBx * VY_PCAT / 2.0 
		+ alphat * alphat * betay * betaz * BZ_PCFT / 2.0 
		- alphaf * alphaf * betay * cf * betaz * sgnBx * sqrt(r) / cs * VZ_MCST / 2.0 
		- alphat * betay * alphaf / sqrt(r) / cs * P_MCST / 2.0 
		+ alphat * betay * alphaf / sqrt(r) / cs * P_MCFT / 2.0 
		+ alphaf * betay * alphat * ct * sqrt(r) / cs * VX_PCST / 2.0 
		+ alphaf * alphaf * betay * betay * cf * sgnBx * sqrt(r) / cs * VY_PCST / 2.0 
		- alphat * alphat * betay * betay * ct * sgnBx * sqrt(r) / cs * VY_MCFT / 2.0 
		- alphat * alphat * betay * ct * betaz * sgnBx * sqrt(r) / cs * VZ_MCFT / 2.0 
		+ alphaf * alphaf * betay * betay * BY_PCST / 2.0 
		- betaz * betay * BZ_MCAT / 2.0 
		+ alphaf * alphaf * betay * betay * BY_MCST / 2.0 
		- betaz * betay * BZ_PCAT / 2.0 
		+ alphat * alphat * betay * betay * BY_PCFT / 2.0 
		+ alphat * betay * alphaf * cf * sqrt(r) / cs * VX_MCFT / 2.0 
		- alphaf * alphaf * betay * betay * cf * sgnBx * sqrt(r) / cs * VY_MCST / 2.0 
		- betaz * sqrt(r) * betay * sgnBx * VZ_PCAT / 2.0 
		+ alphat * alphat * betay * betay * ct * sgnBx * sqrt(r) / cs * VY_PCFT / 2.0 
		+ alphat * alphat * betay * ct * betaz * sgnBx * sqrt(r) / cs * VZ_PCFT / 2.0 
		+ alphat * betay * alphaf * P_PCFT / sqrt(r) / cs / 2.0;

	fbz = alphaf * alphaf * betaz * betay * BY_PCST / 2.0 
		- alphat * alphat * betaz * ct * betay * sgnBx * sqrt(r) / cs * VY_MCFT / 2.0 
		+ betay * betay * BZ_PCAT / 2.0 
		- alphat * alphat * betaz * betaz * ct * sgnBx * sqrt(r) / cs * VZ_MCFT / 2.0 
		+ betay * betay * BZ_MCAT / 2.0 
		+ alphat * alphat * betaz * betay * BY_PCFT / 2.0 
		+ alphat * alphat * betaz * ct * betay * sgnBx * sqrt(r) / cs * VY_PCFT / 2.0 
		- betay * betaz * BY_PCAT / 2.0 
		+ alphaf * alphaf * betaz * betay * BY_MCST / 2.0 
		- alphaf * betaz * alphat * ct * sqrt(r) / cs * VX_MCST / 2.0 
		- betay * betaz * BY_MCAT / 2.0 
		- betay * betay * sqrt(r) * sgnBx * VZ_MCAT / 2.0 
		+ alphat * alphat * betaz * betay * BY_MCFT / 2.0 
		+ alphat * betaz * alphaf * P_PCFT / sqrt(r) / cs / 2.0 
		+ alphat * alphat * betaz * betaz * BZ_PCFT / 2.0 
		+ alphat * alphat * betaz * betaz * BZ_MCFT / 2.0 
		- alphat * betaz * alphaf / sqrt(r) / cs * P_MCST / 2.0 
		+ alphaf * alphaf * betaz * betaz * BZ_PCST / 2.0 
		- alphat * betaz * alphaf / sqrt(r) / cs * P_PCST / 2.0 
		+ alphaf * alphaf * betaz * betaz * cf * sgnBx * sqrt(r) / cs * VZ_PCST / 2.0 
		- alphaf * alphaf * betaz * cf * betay * sgnBx * sqrt(r) / cs * VY_MCST / 2.0 
		+ betay * betay * sqrt(r) * sgnBx * VZ_PCAT / 2.0 
		+ alphaf * alphaf * betaz * betaz * BZ_MCST / 2.0 
		+ alphat * betaz * alphaf / sqrt(r) / cs * P_MCFT / 2.0 
		+ alphat * betaz * alphaf * cf * sqrt(r) / cs * VX_MCFT / 2.0 
		+ alphaf * alphaf * betaz * cf * betay * sgnBx * sqrt(r) / cs * VY_PCST / 2.0 
		+ alphaf * betaz * alphat * ct * sqrt(r) / cs * VX_PCST / 2.0 
		- betay * sqrt(r) * betaz * sgnBx * VY_PCAT / 2.0 
		+ alphat * alphat * betaz * betaz * ct * sgnBx * sqrt(r) / cs * VZ_PCFT / 2.0 
		+ betay * sqrt(r) * betaz * sgnBx * VY_MCAT / 2.0 
		- alphat * betaz * alphaf * cf * sqrt(r) / cs * VX_PCFT / 2.0 
		- alphaf * alphaf * betaz * betaz * cf * sgnBx * sqrt(r) / cs * VZ_MCST / 2.0;

	fp =  r * alphaf * alphaf * cf * VX_MCFT / 2.0 
		- r * alphaf * alphat * ct * betay * sgnBx * VY_MCFT / 2.0 
		- r * alphaf * alphat * ct * betaz * sgnBx * VZ_MCFT / 2.0 
		+ sqrt(r) * alphaf * alphat * betay * BY_MCFT * cs / 2.0 
		+ sqrt(r) * alphaf * alphat * betaz * BZ_MCFT * cs / 2.0 
		+ alphaf * alphaf * P_MCFT / 2.0 
		+ r * alphaf * alphat * ct * betay * sgnBx * VY_PCFT / 2.0 
		+ r * alphaf * alphat * ct * betaz * sgnBx * VZ_PCFT / 2.0 
		+ sqrt(r) * alphaf * alphat * betay * BY_PCFT * cs / 2.0 
		+ sqrt(r) * alphaf * alphat * betaz * BZ_PCFT * cs / 2.0 
		+ alphaf * alphaf * P_PCFT / 2.0 
		+ r * alphat * alphaf * cf * betay * sgnBx * VY_MCST / 2.0 
		+ r * alphat * alphaf * cf * betaz * sgnBx * VZ_MCST / 2.0 
		- sqrt(r) * alphat * alphaf * betay * BY_MCST * cs / 2.0 
		- sqrt(r) * alphat * alphaf * betaz * BZ_MCST * cs / 2.0 
		+ alphat * alphat * P_MCST / 2.0 
		- r * alphat * alphat * ct * VX_PCST / 2.0 
		- r * alphat * alphaf * cf * betay * sgnBx * VY_PCST / 2.0 
		- r * alphat * alphaf * cf * betaz * sgnBx * VZ_PCST / 2.0 
		- sqrt(r) * alphat * alphaf * betay * BY_PCST * cs / 2.0 
		- sqrt(r) * alphat * alphaf * betaz * BZ_PCST * cs / 2.0 
		+ alphat * alphat * P_PCST / 2.0 
		- r * alphaf * alphaf * cf * VX_PCFT / 2.0 
		+ r * alphat * alphat * ct * VX_MCST / 2.0;

	fbx = bx;

}

/*******************************************************************/
#endif