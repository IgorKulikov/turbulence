#ifndef _advec_hpp_
#define _advec_hpp_
/*******************************************************************/
#include "global.hpp"
#include "boundary.hpp"
#include "system.hpp"
#include "parabola.hpp"

real median8(real vcntr, real vpz, real vmz, real vzp, real vzm,
			real vpp, real vpm, real vmp, real vmm)
{
	return 0.0625 * (4.0*vcntr + 2.0*(vpz + vmz + vzp + vzm) + vpp + vpm + vmp + vmm);
}

void advec(mesh *a, real *afin, real *rho, real *vx, real *vy, real *vz)
{
	integer i, k, l;
	real sqrxyz, sqrpx, sqrmx, sqrpy, sqrmy, sqrpz, sqrmz;
	real vxcn, vycn, vzcn, vxpx, vxmx, vypy, vymy, vzpz, vzmz;
	real vpx, vmx, vpy, vmy, vpz, vmz, div;
	real flapx, flamx, flapy, flamy, flapz, flamz;
	
	for(i=1;i<NX-1;i++)
		for(k=1;k<NY-1;k++)
			for(l=1;l<NZ-1;l++)
			{
				// density
				sqrxyz = sqrt(rho[index(i,k,l)]); 
				sqrpx  = sqrt(rho[index(i+1,k,l)]);
				sqrmx  = sqrt(rho[index(i-1,k,l)]);
				sqrpy  = sqrt(rho[index(i,k+1,l)]);
				sqrmy  = sqrt(rho[index(i,k-1,l)]);
				sqrpz  = sqrt(rho[index(i,k,l+1)]);
				sqrmz  = sqrt(rho[index(i,k,l-1)]);

				// velocity
				vxcn = median8(vx[index(i,k,l)],
					vx[index(i,k+1,l)],vx[index(i,k-1,l)],
					vx[index(i,k,l+1)],vx[index(i,k,l-1)],
					vx[index(i,k+1,l+1)],vx[index(i,k+1,l-1)],
					vx[index(i,k-1,l+1)],vx[index(i,k-1,l-1)]);
				vycn = median8(vy[index(i,k,l)],
					vy[index(i+1,k,l)],vy[index(i-1,k,l)],
					vy[index(i,k,l+1)],vy[index(i,k,l-1)],
					vy[index(i+1,k,l+1)],vy[index(i+1,k,l-1)],
					vy[index(i-1,k,l+1)],vy[index(i-1,k,l-1)]);
				vzcn = median8(vz[index(i,k,l)],
					vz[index(i+1,k,l)],vz[index(i-1,k,l)],
					vz[index(i,k+1,l)],vz[index(i,k-1,l)],
					vz[index(i+1,k+1,l)],vz[index(i+1,k-1,l)],
					vz[index(i-1,k+1,l)],vz[index(i-1,k-1,l)]);
				vxpx = median8(vx[index(i+1,k,l)],
					vx[index(i+1,k+1,l)],vx[index(i+1,k-1,l)],
					vx[index(i+1,k,l+1)],vx[index(i+1,k,l-1)],
					vx[index(i+1,k+1,l+1)],vx[index(i+1,k+1,l-1)],
					vx[index(i+1,k-1,l+1)],vx[index(i+1,k-1,l-1)]);
				vxmx = median8(vx[index(i-1,k,l)],
					vx[index(i-1,k+1,l)],vx[index(i-1,k-1,l)],
					vx[index(i-1,k,l+1)],vx[index(i-1,k,l-1)],
					vx[index(i-1,k+1,l+1)],vx[index(i-1,k+1,l-1)],
					vx[index(i-1,k-1,l+1)],vx[index(i-1,k-1,l-1)]);
				vypy = median8(vy[index(i,k+1,l)],
					vy[index(i+1,k+1,l)],vy[index(i-1,k+1,l)],
					vy[index(i,k+1,l+1)],vy[index(i,k+1,l-1)],
					vy[index(i+1,k+1,l+1)],vy[index(i+1,k+1,l-1)],
					vy[index(i-1,k+1,l+1)],vy[index(i-1,k+1,l-1)]);
				vymy = median8(vy[index(i,k-1,l)],
					vy[index(i+1,k-1,l)],vy[index(i-1,k-1,l)],
					vy[index(i,k-1,l+1)],vy[index(i,k-1,l-1)],
					vy[index(i+1,k-1,l+1)],vy[index(i+1,k-1,l-1)],
					vy[index(i-1,k-1,l+1)],vy[index(i-1,k-1,l-1)]);
				vzpz = median8(vz[index(i,k,l+1)],
					vz[index(i+1,k,l+1)],vz[index(i-1,k,l+1)],
					vz[index(i,k+1,l+1)],vz[index(i,k-1,l+1)],
					vz[index(i+1,k+1,l+1)],vz[index(i+1,k-1,l+1)],
					vz[index(i-1,k+1,l+1)],vz[index(i-1,k-1,l+1)]);
				vzmz = median8(vz[index(i,k,l-1)],
					vz[index(i+1,k,l-1)],vz[index(i-1,k,l-1)],
					vz[index(i,k+1,l-1)],vz[index(i,k-1,l-1)],
					vz[index(i+1,k+1,l-1)],vz[index(i+1,k-1,l-1)],
					vz[index(i-1,k+1,l-1)],vz[index(i-1,k-1,l-1)]);

				// characteristic
				vpx = (sqrpx * vxpx + sqrxyz * vxcn)/(sqrxyz + sqrpx);
				vmx = (sqrmx * vxmx + sqrxyz * vxcn)/(sqrxyz + sqrmx);
				vpy = (sqrpy * vypy + sqrxyz * vycn)/(sqrxyz + sqrpy);
				vmy = (sqrmy * vymy + sqrxyz * vycn)/(sqrxyz + sqrmy);
				vpz = (sqrpz * vzpz + sqrxyz * vzcn)/(sqrxyz + sqrpz);
				vmz = (sqrmz * vzmz + sqrxyz * vzcn)/(sqrxyz + sqrmz);

				// Computational flux
				if(vpx > 0) 
					flapx = vpx * a[index(i,k,l)].parx.right_parabola(vpx*tau,h);
				else
					flapx = vpx * a[index(i+1,k,l)].parx.left_parabola(-vpx*tau,h);
							
				if(vmx > 0)
					flamx = vmx * a[index(i-1,k,l)].parx.right_parabola(vmx*tau,h);
				else
					flamx = vmx * a[index(i,k,l)].parx.left_parabola(-vmx*tau,h);

				if(vpy > 0)
					flapy = vpy * a[index(i,k,l)].pary.right_parabola(vpy*tau,h);
				else
					flapy = vpy * a[index(i,k+1,l)].pary.left_parabola(-vpy*tau,h);
				
				if(vmy > 0)
					flamy = vmy * a[index(i,k-1,l)].pary.right_parabola(vmy*tau,h);
				else
					flamy = vmy * a[index(i,k,l)].pary.left_parabola(-vmy*tau,h);

				if(vpz > 0)
					flapz = vpz * a[index(i,k,l)].parz.right_parabola(vpz*tau,h);
				else
					flapz = vpz * a[index(i,k,l+1)].parz.left_parabola(-vpz*tau,h);

				if(vmz > 0)
					flamz = vmz * a[index(i,k,l-1)].parz.right_parabola(vmz*tau,h);
				else
					flamz = vmz * a[index(i,k,l)].parz.left_parabola(-vmz*tau,h);

				// Equation
				div = (flapx - flamx)/h + (flapy - flamy)/h + (flapz - flamz)/h;
				afin[index(i,k,l)] = a[index(i,k,l)].u - tau * div;
			}

	boundary(afin);
}

/*******************************************************************/
#endif