/*  AstroPhi 2.0  */
#include "init.hpp"
#include "boundary.hpp"
#include "build_parabola.hpp"
#include "parallel.hpp"
#include "global.hpp"
#include "parabola.hpp"
#include "mesh.hpp"
#include "system.hpp"
#include "problem.hpp"
#include "save.hpp"
#include "cfl.hpp"
#include "imp2vel.hpp"
#include "euler.hpp"
#include "lagrange.hpp"
#include "poisson.hpp"
#include "entropy.hpp"
#include "analysis.hpp"

int main(int argc, char *argv [])
{	
	// file with conservationa laws
	FILE* fout = NULL;
	
	// Current time
	real timer = 0.0, tend;

	// Start MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0) fout = fopen("const.log","w");

	// Initial computational
	parallel_init();
	hydro_init();

	// Load problem
	load_problem();

	// Save init
	save(R,Bx,By,timer);
	
	for(integer i=0 ; i<iTimeCount ; i++)
	{
		tend = print_tau * (i+1);
		while(timer < tend)
		{
			// MHD
			tau = computational_tau(timer,tend);
			timer += tau;
			laplass();
			eulerian_stage();
			imp2vel(R_Next,RUx_Next,RUy_Next,RUz_Next,Vx,Vy,Vz);
			lagrange_stage();
			imp2vel(R,RUx,RUy,RUz,Vx,Vy,Vz);
			
			// Entropy correction
			entropy(Vx,Vy,Vz,Bx,By,Bz,P,RE,R);
			
			save_log(fout,timer);
		}
		save(R,Bx,By,timer);
		if(rank == 0) printf("Time = %lf\n",timer);
	}
		analysis(R,Vx,Vy,Vz,Bx,By,Bz,timer);
	
	// Destroy
	if(rank == 0) fclose(fout);
	hydro_destroy();
	parallel_destroy();
	MPI_Finalize();
	return 0;	
}
