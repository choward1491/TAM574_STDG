#include "simulation.h"

void Simulation::run(int numx,int numt, string fileu, string fileq,string filex)
{
	if( !grid.empty() )
	{
		int rows = grid.Rows(), cols = grid.Cols();
		
		for(int i = 0; i < rows; i++)
		{
			cout<<"The simulation is "<<(double)i*100/(double)rows<<"% Done Computing\n";
			grid.SolveConstants(i); 
		}
		//cout<<"The simulation is 100% Done Computing\n";

		//cout<<"P = "<<grid(1,1).getP()<<", (Nx,Nt) = ("<<grid.Cols()<<", "<<grid.Rows()<<"), zeta_u = "<<grid.convergence_u()<<" & zeta_q = "<<grid.convergence_q()<<endl;

		//grid.writeSolution2File(numx,numt,fileu,fileq,filex);
		//grid.writeEndSolution2File(numx,numt,fileu,fileq,filex);
		//cout<<"Done writing the solution data to file\n";
	}
}