#ifndef SIMULATION_H
#define SIMULATION_H

#include "cell.h"
#include "grid.hpp"
#include "postprocessing.h"

class Simulation
{
	public:
		Simulation(int rows, int cols, 
			int poly_order, 
			double width, double height,
			double (*iu)(double x),double (*iq)(double x),
			Vector constants)
		{
			grid.change(rows,cols,poly_order,width,height,constants);
			grid.initU = iu;
			grid.initQ = iq;
		}

    void run(int numx,int numt,
             string full_fileu, string full_fileq,string full_filex,
             string final_fileu, string final_fileq,string final_filex);
	private:
		Grid grid;
};


#endif
