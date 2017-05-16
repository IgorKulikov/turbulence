#ifndef _euler_hpp_
#define _euler_hpp_
/*******************************************************************/
#include "global.hpp"
#include "boundary.hpp"
#include "riemann.hpp"
#include "build_parabola.hpp"

void eulerian_stage()
{
	integer i, k, l;
	
	real fxvxp, fxvyp, fxvzp, fxpp, fxbxp, fxbyp, fxbzp;
	real fxvxm, fxvym, fxvzm, fxpm, fxbxm, fxbym, fxbzm;
	real fyvyp, fyvzp, fyvxp, fypp, fybyp, fybzp, fybxp;
	real fyvym, fyvzm, fyvxm, fypm, fybym, fybzm, fybxm;
	real fzvzp, fzvxp, fzvyp, fzpp, fzbzp, fzbxp, fzbyp;
	real fzvzm, fzvxm, fzvym, fzpm, fzbzm, fzbxm, fzbym;

	real tpresxp, tpresxm, tpresyp, tpresym, tpreszp, tpreszm;
	real dotxp, dotxm, dotyp, dotym, dotzp, dotzm;
	real gradfix, gradfiy, gradfiz, divv, divpv;
	
	build_parabola(Vx,parvx);
	build_parabola(Vy,parvy);
	build_parabola(Vz,parvz);
	build_parabola(P,parp);
	build_parabola(Bx,parbx);
	build_parabola(By,parby);
	build_parabola(Bz,parbz);

	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				// Riemann solver for Hydrodynamics equations
				mhd_riemann(	/* (i,k,l)-(i+1,k,l) interface */
					R[index(i+1,k,l)],			R[index(i,k,l)], 
					P[index(i+1,k,l)],			P[index(i,k,l)], 
					Vx[index(i+1,k,l)],			Vx[index(i,k,l)], 
					Vy[index(i+1,k,l)],			Vy[index(i,k,l)], 
					Vz[index(i+1,k,l)],			Vz[index(i,k,l)],
					Bx[index(i+1,k,l)],			Bx[index(i,k,l)], 
					By[index(i+1,k,l)],			By[index(i,k,l)], 
					Bz[index(i+1,k,l)],			Bz[index(i,k,l)],
					parp[index(i+1,k,l)].parx,	parp[index(i,k,l)].parx,
					parvx[index(i+1,k,l)].parx,	parvx[index(i,k,l)].parx,
					parvy[index(i+1,k,l)].parx,	parvy[index(i,k,l)].parx,
					parvz[index(i+1,k,l)].parx,	parvz[index(i,k,l)].parx,
					parby[index(i+1,k,l)].parx,	parby[index(i,k,l)].parx,
					parbz[index(i+1,k,l)].parx,	parbz[index(i,k,l)].parx,
					fxvxp, fxvyp, fxvzp, fxpp, fxbxp, fxbyp, fxbzp);

				mhd_riemann(	/* (i-1,k,l)-(i,k,l) interface */
					R[index(i,k,l)],			R[index(i-1,k,l)], 
					P[index(i,k,l)],			P[index(i-1,k,l)], 
					Vx[index(i,k,l)],			Vx[index(i-1,k,l)], 
					Vy[index(i,k,l)],			Vy[index(i-1,k,l)], 
					Vz[index(i,k,l)],			Vz[index(i-1,k,l)],
					Bx[index(i,k,l)],			Bx[index(i-1,k,l)], 
					By[index(i,k,l)],			By[index(i-1,k,l)], 
					Bz[index(i,k,l)],			Bz[index(i-1,k,l)],
					parp[index(i,k,l)].parx,	parp[index(i-1,k,l)].parx,
					parvx[index(i,k,l)].parx,	parvx[index(i-1,k,l)].parx,
					parvy[index(i,k,l)].parx,	parvy[index(i-1,k,l)].parx,
					parvz[index(i,k,l)].parx,	parvz[index(i-1,k,l)].parx,
					parby[index(i,k,l)].parx,	parby[index(i-1,k,l)].parx,
					parbz[index(i,k,l)].parx,	parbz[index(i-1,k,l)].parx,
					fxvxm, fxvym, fxvzm, fxpm, fxbxm, fxbym, fxbzm);

				mhd_riemann(	/* (i,k,l)-(i,k+1,l) interface */
					R[index(i,k+1,l)],			R[index(i,k,l)], 
					P[index(i,k+1,l)],			P[index(i,k,l)], 
					Vy[index(i,k+1,l)],			Vy[index(i,k,l)], 
					Vz[index(i,k+1,l)],			Vz[index(i,k,l)], 
					Vx[index(i,k+1,l)],			Vx[index(i,k,l)],
					By[index(i,k+1,l)],			By[index(i,k,l)], 
					Bz[index(i,k+1,l)],			Bz[index(i,k,l)], 
					Bx[index(i,k+1,l)],			Bx[index(i,k,l)],
					parp[index(i,k+1,l)].pary,	parp[index(i,k,l)].pary,
					parvy[index(i,k+1,l)].pary,	parvy[index(i,k,l)].pary,
					parvz[index(i,k+1,l)].pary,	parvz[index(i,k,l)].pary,
					parvx[index(i,k+1,l)].pary,	parvx[index(i,k,l)].pary,
					parbz[index(i,k+1,l)].pary,	parbz[index(i,k,l)].pary,
					parbx[index(i,k+1,l)].pary,	parbx[index(i,k,l)].pary,
					fyvyp, fyvzp, fyvxp, fypp, fybyp, fybzp, fybxp);

				mhd_riemann(	/* (i,k-1,l)-(i,k,l) interface */
					R[index(i,k,l)],			R[index(i,k-1,l)], 
					P[index(i,k,l)],			P[index(i,k-1,l)], 
					Vy[index(i,k,l)],			Vy[index(i,k-1,l)], 
					Vz[index(i,k,l)],			Vz[index(i,k-1,l)], 
					Vx[index(i,k,l)],			Vx[index(i,k-1,l)],
					By[index(i,k,l)],			By[index(i,k-1,l)], 
					Bz[index(i,k,l)],			Bz[index(i,k-1,l)], 
					Bx[index(i,k,l)],			Bx[index(i,k-1,l)],
					parp[index(i,k,l)].pary,	parp[index(i,k-1,l)].pary,
					parvy[index(i,k,l)].pary,	parvy[index(i,k-1,l)].pary,
					parvz[index(i,k,l)].pary,	parvz[index(i,k-1,l)].pary,
					parvx[index(i,k,l)].pary,	parvx[index(i,k-1,l)].pary,
					parbz[index(i,k,l)].pary,	parbz[index(i,k-1,l)].pary,
					parbx[index(i,k,l)].pary,	parbx[index(i,k-1,l)].pary,
					fyvym, fyvzm, fyvxm, fypm, fybym, fybzm, fybxm);

				mhd_riemann(	/* (i,k,l)-(i,k,l+1) interface */
					R[index(i,k,l+1)],			R[index(i,k,l)], 
					P[index(i,k,l+1)],			P[index(i,k,l)], 
					Vz[index(i,k,l+1)],			Vz[index(i,k,l)], 
					Vx[index(i,k,l+1)],			Vx[index(i,k,l)], 
					Vy[index(i,k,l+1)],			Vy[index(i,k,l)],
					Bz[index(i,k,l+1)],			Bz[index(i,k,l)], 
					Bx[index(i,k,l+1)],			Bx[index(i,k,l)], 
					By[index(i,k,l+1)],			By[index(i,k,l)],
					parp[index(i,k,l+1)].parz,	parp[index(i,k,l)].parz,
					parvz[index(i,k,l+1)].parz,	parvz[index(i,k,l)].parz,
					parvx[index(i,k,l+1)].parz,	parvx[index(i,k,l)].parz,
					parvy[index(i,k,l+1)].parz,	parvy[index(i,k,l)].parz,
					parbx[index(i,k,l+1)].parz,	parbx[index(i,k,l)].parz,
					parby[index(i,k,l+1)].parz,	parby[index(i,k,l)].parz,
					fzvzp, fzvxp, fzvyp, fzpp, fzbzp, fzbxp, fzbyp);

				mhd_riemann(	/* (i,k,l-1)-(i,k,l) interface */
					R[index(i,k,l)],			R[index(i,k,l-1)], 
					P[index(i,k,l)],			P[index(i,k,l-1)], 
					Vz[index(i,k,l)],			Vz[index(i,k,l-1)], 
					Vx[index(i,k,l)],			Vx[index(i,k,l-1)], 
					Vy[index(i,k,l)],			Vy[index(i,k,l-1)],
					Bz[index(i,k,l)],			Bz[index(i,k,l-1)], 
					Bx[index(i,k,l)],			Bx[index(i,k,l-1)], 
					By[index(i,k,l)],			By[index(i,k,l-1)],
					parp[index(i,k,l)].parz,	parp[index(i,k,l-1)].parz,
					parvz[index(i,k,l)].parz,	parvz[index(i,k,l-1)].parz,
					parvx[index(i,k,l)].parz,	parvx[index(i,k,l-1)].parz,
					parvy[index(i,k,l)].parz,	parvy[index(i,k,l-1)].parz,
					parbx[index(i,k,l)].parz,	parbx[index(i,k,l-1)].parz,
					parby[index(i,k,l)].parz,	parby[index(i,k,l-1)].parz,
					fzvzm, fzvxm, fzvym, fzpm, fzbzm, fzbxm, fzbym);

				// Equation
				gradfix = (Fi[index(i+1,k,l)] - Fi[index(i-1,k,l)])/(2.0*h);
				gradfiy = (Fi[index(i,k+1,l)] - Fi[index(i,k-1,l)])/(2.0*h);
				gradfiz = (Fi[index(i,k,l+1)] - Fi[index(i,k,l-1)])/(2.0*h);

				tpresxp = fxpp + 0.5 * (fxbxp*fxbxp + fxbyp*fxbyp + fxbzp*fxbzp);
				tpresxm = fxpm + 0.5 * (fxbxm*fxbxm + fxbym*fxbym + fxbzm*fxbzm);
				tpresyp = fypp + 0.5 * (fybxp*fybxp + fybyp*fybyp + fybzp*fybzp);
				tpresym = fypm + 0.5 * (fybxm*fybxm + fybym*fybym + fybzm*fybzm);
				tpreszp = fzpp + 0.5 * (fzbxp*fzbxp + fzbyp*fzbyp + fzbzp*fzbzp);
				tpreszm = fzpm + 0.5 * (fzbxm*fzbxm + fzbym*fzbym + fzbzm*fzbzm);

				dotxp = fxbxp*fxvxp + fxbyp*fxvyp + fxbzp*fxvzp;
				dotxm = fxbxm*fxvxm + fxbym*fxvym + fxbzm*fxvzm;
				dotyp = fybxp*fyvxp + fybyp*fyvyp + fybzp*fyvzp;
				dotym = fybxm*fyvxm + fybym*fyvym + fybzm*fyvzm;
				dotzp = fzbxp*fzvxp + fzbyp*fzvyp + fzbzp*fzvzp;
				dotzm = fzbxm*fzvxm + fzbym*fzvym + fzbzm*fzvzm;

				divv  = (fxvxp - fxvxm)/h + (fyvyp - fyvym)/h + (fzvzp - fzvzm)/h;
				divpv = (fxvxp*tpresxp - fxvxm*tpresxm)/h + 
						(fyvyp*tpresyp - fyvym*tpresym)/h + 
						(fzvzp*tpreszp - fzvzm*tpreszm)/h;

				R_Next[index(i,k,l)] = R[index(i,k,l)];
				
				RUx_Next[index(i,k,l)] = RUx[index(i,k,l)] 
					- tau * (tpresxp - tpresxm)/h
					+ tau * (fxbxp*fxbxp - fxbxm*fxbxm)/h
					+ tau * (fybxp*fybyp - fybxm*fybym)/h
					+ tau * (fzbxp*fzbzp - fzbxm*fzbzm)/h
					- tau * R[index(i,k,l)] * gradfix;

				RUy_Next[index(i,k,l)] = RUy[index(i,k,l)] 
					- tau * (tpresyp - tpresym)/h
					+ tau * (fxbyp*fxbxp - fxbym*fxbxm)/h
					+ tau * (fybyp*fybyp - fybym*fybym)/h
					+ tau * (fzbyp*fzbzp - fzbym*fzbzm)/h
					- tau * R[index(i,k,l)] * gradfiy;

				RUz_Next[index(i,k,l)] = RUz[index(i,k,l)] 
					- tau * (tpreszp - tpreszm)/h
					+ tau * (fxbzp*fxbxp - fxbzm*fxbxm)/h
					+ tau * (fybzp*fybyp - fybzm*fybym)/h
					+ tau * (fzbzp*fzbzp - fzbzm*fzbzm)/h
					- tau * R[index(i,k,l)] * gradfiz;

				P_Next[index(i,k,l)] = P[index(i,k,l)] 
					- tau * (GAMMA - 1.0) * P[index(i,k,l)] * divv;

				RE_Next[index(i,k,l)] = RE[index(i,k,l)] 
					- tau * divpv
					+ tau * (fxbxp*dotxp - fxbxm*dotxm)/h
					+ tau * (fybyp*dotyp - fybym*dotym)/h
					+ tau * (fzbzp*dotzp - fzbzm*dotzm)/h
					- RUx[index(i,k,l)] * gradfix
					- RUy[index(i,k,l)] * gradfiy
					- RUz[index(i,k,l)] * gradfiz;

				Bx_Next[index(i,k,l)] = Bx[index(i,k,l)] 
					- tau * ((fybxp*fyvyp-fybyp*fyvxp) - (fybxm*fyvym-fybym*fyvxm))/h
					- tau * ((fzbxp*fzvzp-fzbzp*fzvxp) - (fzbxm*fzvzm-fzbzm*fzvxm))/h;
				By_Next[index(i,k,l)] = By[index(i,k,l)] 
					- tau * ((fxbyp*fxvxp-fxbxp*fxvyp) - (fxbym*fxvxm-fxbxm*fxvym))/h
					- tau * ((fzbyp*fzvzp-fzbzp*fzvyp) - (fzbym*fzvzm-fzbzm*fzvym))/h;
				Bz_Next[index(i,k,l)] = Bz[index(i,k,l)]
					- tau * ((fxbzp*fxvxp-fxbxp*fxvzp) - (fxbzm*fxvxm-fxbxm*fxvzm))/h
					- tau * ((fybzp*fyvyp-fybyp*fyvzp) - (fybzm*fyvym-fybym*fyvzm))/h;
			}

	boundary(R_Next);
	boundary(RUx_Next);
	boundary(RUy_Next);
	boundary(RUz_Next);
	boundary(P_Next);
	boundary(RE_Next);
	boundary(Bx_Next);
	boundary(By_Next);
	boundary(Bz_Next);	
}

/*******************************************************************/
#endif