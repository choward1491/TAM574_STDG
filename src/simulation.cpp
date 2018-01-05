#include "simulation.h"

void Simulation::run(int numx,int numt,
                     string full_fileu, string full_fileq,string full_filex,
                     string final_fileu, string final_fileq,string final_filex)
{
	if( !grid.empty() )
	{
		int rows = grid.Rows(), cols = grid.Cols();
		
		for(int i = 0; i < rows; i++)
		{
			cout<<"The simulation is "<<(double)i*100/(double)rows<<"% Done Computing\n";
			grid.SolveConstants(i); 
		}
		cout<<"The simulation is 100% Done Computing\n";
        cout << "Compute convergence related factors.." << endl;
		cout<<"P = "<<grid(1,1).getP()<<", (Nx,Nt) = ("<<grid.Cols()<<", "<<grid.Rows()<<"), zeta_u = "<<grid.convergence_u()<<" & zeta_q = "<<grid.convergence_q()<<endl;
        cout << "Writing data to file.." << endl;
		grid.writeSolution2File(numx,numt,full_fileu,full_fileq,full_filex);
		grid.writeEndSolution2File(numx,numt,final_fileu,final_fileq,final_filex);
		cout<<"Done writing the solution data to file.\n";
	}
}
