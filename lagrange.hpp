#ifndef _lagrange_hpp_
#define _lagrange_hpp_
/*******************************************************************/
#include "global.hpp"
#include "advec.hpp"

void lagrange_stage()
{
	integer i, k, l;
	
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				Bx[index(i,k,l)] = Bx_Next[index(i,k,l)];
				By[index(i,k,l)] = By_Next[index(i,k,l)];
				Bz[index(i,k,l)] = Bz_Next[index(i,k,l)];
			}
	boundary(Bx);
	boundary(By);
	boundary(Bz);	


	build_parabola(RUx_Next,parp);
	advec(parp,RUx,R_Next,Vx,Vy,Vz);

	build_parabola(RUy_Next,parp);
	advec(parp,RUy,R_Next,Vx,Vy,Vz);

	build_parabola(RUz_Next,parp);
	advec(parp,RUz,R_Next,Vx,Vy,Vz);

	build_parabola(P_Next,parp);
	advec(parp,P,R_Next,Vx,Vy,Vz);

	build_parabola(RE_Next,parp);
	advec(parp,RE,R_Next,Vx,Vy,Vz);

	build_parabola(R_Next,parp);
	advec(parp,R,R_Next,Vx,Vy,Vz);

}

/*******************************************************************/
#endif