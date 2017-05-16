#ifndef _boundary_hpp_
#define _boundary_hpp_
/*******************************************************************/
#include "mesh.hpp"
#include "system.hpp"
#include "parallel.hpp"
#include "global.hpp"

/* boundary conditions for mesh */
void micro_move(mesh *a, real *buffer, integer i, integer direct)
{
	integer k, l;
	for(k=0;k<NY;k++)
		for(l=0;l<NZ;l++)
			if(direct == 1)
			{
				buffer[k*16*NZ+16*l+0]  = a[index(i,k,l)].u;
				buffer[k*16*NZ+16*l+1]  = a[index(i,k,l)].parx.ddq;
				buffer[k*16*NZ+16*l+2]  = a[index(i,k,l)].parx.dq;
				buffer[k*16*NZ+16*l+3]  = a[index(i,k,l)].parx.dq6;
				buffer[k*16*NZ+16*l+4]  = a[index(i,k,l)].parx.ull;
				buffer[k*16*NZ+16*l+5]  = a[index(i,k,l)].parx.ulr;
				buffer[k*16*NZ+16*l+6]  = a[index(i,k,l)].pary.ddq;
				buffer[k*16*NZ+16*l+7]  = a[index(i,k,l)].pary.dq;
				buffer[k*16*NZ+16*l+8]  = a[index(i,k,l)].pary.dq6;
				buffer[k*16*NZ+16*l+9]  = a[index(i,k,l)].pary.ull;
				buffer[k*16*NZ+16*l+10] = a[index(i,k,l)].pary.ulr;
				buffer[k*16*NZ+16*l+11] = a[index(i,k,l)].parz.ddq;
				buffer[k*16*NZ+16*l+12] = a[index(i,k,l)].parz.dq;
				buffer[k*16*NZ+16*l+13] = a[index(i,k,l)].parz.dq6;
				buffer[k*16*NZ+16*l+14] = a[index(i,k,l)].parz.ull;
				buffer[k*16*NZ+16*l+15] = a[index(i,k,l)].parz.ulr;
			}
			else
			{
				a[index(i,k,l)].u        = buffer[k*16*NZ+16*l+0];  
				a[index(i,k,l)].parx.ddq = buffer[k*16*NZ+16*l+1]; 
				a[index(i,k,l)].parx.dq  = buffer[k*16*NZ+16*l+2];  
				a[index(i,k,l)].parx.dq6 = buffer[k*16*NZ+16*l+3];  
				a[index(i,k,l)].parx.ull = buffer[k*16*NZ+16*l+4];  
				a[index(i,k,l)].parx.ulr = buffer[k*16*NZ+16*l+5];  
				a[index(i,k,l)].pary.ddq = buffer[k*16*NZ+16*l+6];  
				a[index(i,k,l)].pary.dq  = buffer[k*16*NZ+16*l+7];  
				a[index(i,k,l)].pary.dq6 = buffer[k*16*NZ+16*l+8];  
				a[index(i,k,l)].pary.ull = buffer[k*16*NZ+16*l+9];  
				a[index(i,k,l)].pary.ulr = buffer[k*16*NZ+16*l+10]; 
				a[index(i,k,l)].parz.ddq = buffer[k*16*NZ+16*l+11]; 
				a[index(i,k,l)].parz.dq  = buffer[k*16*NZ+16*l+12]; 
				a[index(i,k,l)].parz.dq6 = buffer[k*16*NZ+16*l+13]; 
				a[index(i,k,l)].parz.ull = buffer[k*16*NZ+16*l+14]; 
				a[index(i,k,l)].parz.ulr = buffer[k*16*NZ+16*l+15]; 
			}
}

void boundary_mesh(mesh *a)
{
	integer i, k, l;

	// boundary conditions on x axis
	if(size == 1)
	{
		for(k=0;k<NY;k++)
			for(l=0;l<NZ;l++)
			{
				a[index(0,k,l)]    = a[index(NX-2,k,l)];
				a[index(NX-1,k,l)] = a[index(1,k,l)];
			}
	}
	else
	{
	// Communications to right
		if(rank == 0)
		{
			micro_move(a,buffer,NX - 2, 1);
			MPI_Send(buffer,16*NY*NZ,MPI_DOUBLE,(rank+1+size)%size,TAGTORIGHT,MPI_COMM_WORLD);
			MPI_Recv(buffer,16*NY*NZ,MPI_DOUBLE,(rank-1+size)%size,TAGTORIGHT,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			micro_move(a,buffer, 0, 0);
		}
		else
		{
			MPI_Recv(buffer,16*NY*NZ,MPI_DOUBLE,(rank-1+size)%size,TAGTORIGHT,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			micro_move(a,buffer, 0, 0);
			micro_move(a,buffer,NX - 2, 1);
			MPI_Send(buffer,16*NY*NZ,MPI_DOUBLE,(rank+1+size)%size,TAGTORIGHT,MPI_COMM_WORLD);
		}
	// Communications to left
		if(rank == size-1)
		{
			micro_move(a,buffer, 1, 1);
			MPI_Send(buffer,16*NY*NZ,MPI_DOUBLE,(rank-1+size)%size,TAGTOLEFT,MPI_COMM_WORLD);
			MPI_Recv(buffer,16*NY*NZ,MPI_DOUBLE,(rank+1+size)%size,TAGTOLEFT,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			micro_move(a,buffer, NX - 1, 0);
		}
		else
		{
			MPI_Recv(buffer,16*NY*NZ,MPI_DOUBLE,(rank+1+size)%size,TAGTOLEFT,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			micro_move(a,buffer, NX - 1, 0);
			micro_move(a,buffer, 1, 1);
			MPI_Send(buffer,16*NY*NZ,MPI_DOUBLE,(rank-1+size)%size,TAGTOLEFT,MPI_COMM_WORLD);
		}
	}

	// boundary conditions on y axis
	for(i=0;i<NX;i++)
		for(l=0;l<NZ;l++)
			{
				a[index(i,0,l)]    = a[index(i,NY-2,l)];
				a[index(i,NY-1,l)] = a[index(i,1,l)];
			}
	
	// boundary conditions on z axis
	for(i=0;i<NX;i++)
		for(k=0;k<NY;k++)
			{
				a[index(i,k,0)]    = a[index(i,k,NZ-2)];
				a[index(i,k,NZ-1)] = a[index(i,k,1)];
			}

}

/* Boundary condition */
void micro_move_real(real *a, real *buffer, integer i, integer direct)
{
	integer k, l;
	for(k=0;k<NY;k++)
		for(l=0;l<NZ;l++)
		{
			if(direct == 1)
				buffer[k*NZ+l]  = a[index(i,k,l)];
			else
				a[index(i,k,l)] = buffer[k*NZ+l];
		}
}

void boundary(real *a)
{
	integer i, k, l;

	// boundary conditions on x axis
	if(size == 1)
	{
		for(k=0;k<NY;k++)
			for(l=0;l<NZ;l++)
			{
				a[index(0,k,l)]    = a[index(NX-2,k,l)];
				a[index(NX-1,k,l)] = a[index(1,k,l)];
			}
	}
	else
	{
	// Communications to right
		if(rank == 0)
		{
			micro_move_real(a,buffer,NX - 2, 1);
			MPI_Send(buffer,NY*NZ,MPI_DOUBLE,(rank+1+size)%size,TAGTORIGHT,MPI_COMM_WORLD);
			MPI_Recv(buffer,NY*NZ,MPI_DOUBLE,(rank-1+size)%size,TAGTORIGHT,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			micro_move_real(a,buffer, 0, 0);
		}
		else
		{
			MPI_Recv(buffer,NY*NZ,MPI_DOUBLE,(rank-1+size)%size,TAGTORIGHT,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			micro_move_real(a,buffer, 0, 0);
			micro_move_real(a,buffer,NX - 2, 1);
			MPI_Send(buffer,NY*NZ,MPI_DOUBLE,(rank+1+size)%size,TAGTORIGHT,MPI_COMM_WORLD);
		}
	// Communications to left
		if(rank == size-1)
		{
			micro_move_real(a,buffer, 1, 1);
			MPI_Send(buffer,NY*NZ,MPI_DOUBLE,(rank-1+size)%size,TAGTOLEFT,MPI_COMM_WORLD);
			MPI_Recv(buffer,NY*NZ,MPI_DOUBLE,(rank+1+size)%size,TAGTOLEFT,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			micro_move_real(a,buffer, NX - 1, 0);
		}
		else
		{
			MPI_Recv(buffer,NY*NZ,MPI_DOUBLE,(rank+1+size)%size,TAGTOLEFT,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			micro_move_real(a,buffer, NX - 1, 0);
			micro_move_real(a,buffer, 1, 1);
			MPI_Send(buffer,NY*NZ,MPI_DOUBLE,(rank-1+size)%size,TAGTOLEFT,MPI_COMM_WORLD);
		}
	}

	// boundary conditions on y axis
	for(i=0;i<NX;i++)
		for(l=0;l<NZ;l++)
			{
				a[index(i,0,l)]    = a[index(i,NY-2,l)];
				a[index(i,NY-1,l)] = a[index(i,1,l)];
			}
	
	// boundary conditions on z axis
	for(i=0;i<NX;i++)
		for(k=0;k<NY;k++)
			{
				a[index(i,k,0)]    = a[index(i,k,NZ-2)];
				a[index(i,k,NZ-1)] = a[index(i,k,1)];
			}

}

/*******************************************************************/
#endif