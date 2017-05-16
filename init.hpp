#ifndef _init_hpp_
#define _init_hpp_
/*******************************************************************/
#include "parabola.hpp"
#include "system.hpp"
#include "global.hpp"

void hydro_init()
{
	R = new real[NX*NY*NZ];
	R_Next = new real[NX*NY*NZ];
	RUx = new real[NX*NY*NZ];
	RUy = new real[NX*NY*NZ];
	RUz = new real[NX*NY*NZ];
	RUx_Next = new real[NX*NY*NZ];
	RUy_Next = new real[NX*NY*NZ];
	RUz_Next = new real[NX*NY*NZ];
	Vx = new real[NX*NY*NZ];
	Vy = new real[NX*NY*NZ];
	Vz = new real[NX*NY*NZ];
	P = new real[NX*NY*NZ];
	RE = new real[NX*NY*NZ];
	RE_Next = new real[NX*NY*NZ];
	Bx = new real[NX*NY*NZ];
	By = new real[NX*NY*NZ];
	Bz = new real[NX*NY*NZ];
	Bx_Next = new real[NX*NY*NZ];
	By_Next = new real[NX*NY*NZ];
	Bz_Next = new real[NX*NY*NZ];
	Fi = new real[NX*NY*NZ];
	P_Next = new real[NX*NY*NZ];
	parvx = new mesh[NX*NY*NZ];
	parvy = new mesh[NX*NY*NZ];
	parvz = new mesh[NX*NY*NZ];
	parp = new mesh[NX*NY*NZ];
	parbx = new mesh[NX*NY*NZ];
	parby = new mesh[NX*NY*NZ];
	parbz = new mesh[NX*NY*NZ];
}

void hydro_destroy()
{
	delete R;
	delete R_Next;
	delete RUx;
	delete RUy;
	delete RUz;
	delete RUx_Next;
	delete RUy_Next;
	delete RUz_Next;	 
	delete Vx;
	delete Vy;
	delete Vz;
	delete P;
	delete RE;
	delete RE_Next;
	delete Bx;
	delete By;
	delete Bz;
	delete Bx_Next;
	delete By_Next;
	delete Bz_Next;
	delete P_Next;
	delete parvx;
	delete parvy;
	delete parvz;
	delete parp;
	delete parbx;
	delete parby;
	delete parbz;
}

/*******************************************************************/
#endif