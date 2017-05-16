#ifndef _mesh_hpp_
#define _mesh_hpp_
/*******************************************************************/
#include "parabola.hpp"
#include "system.hpp"

// type of mesh
struct mesh
{
	// value on mesh
	real u;
	// parabols
	parabola parx, pary, parz;
};

/*******************************************************************/
#endif