#include "vector.h"
#include "matrix.h"
#include "geometry.h"
#include "cell.h"
#include "calculus.h"
#include <vector>
#include "shapefunctions.h"
#include "simulation.h"

#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214

double FunctionU(double x)
{
	return sin(x);
}

double FunctionQ(double x)
{
	return -cos(x);
}

int main( void)
{
    
	Vector consts(4); 
	double T = 5, L = 2*PI;
    double C = 1.0/256.0, a = 2*PI, c = 2*PI;
    double tau_c = 1.0/(4*PI*PI), tau = C*tau_c;
    int p  = 3;
	int Nx = 100, Nt = 10*T/(PI*PI*tau);
    consts[0] = a;      //a
    consts[1] = c;      //c
    consts[2] = tau;    //tau
    consts[3] = 1;      //kappa

    // construct simulations
	Simulation test(Nt,Nx,p,L,T,&FunctionU,&FunctionQ,consts);
    
    // init output solution files
    char filename[256] = {'\0'};
    sprintf(filename, "u_p%i.txt",p);       string ufile = std::string(filename);
    sprintf(filename, "q_p%i.txt",p);       string qfile = std::string(filename);
    sprintf(filename, "x_p%i.txt",p);       string xfile = std::string(filename);
    sprintf(filename, "final_u_p%i.txt",p); string fufile = std::string(filename);
    sprintf(filename, "final_q_p%i.txt",p); string fqfile = std::string(filename);
    sprintf(filename, "final_x_p%i.txt",p); string fxfile = std::string(filename);
    
    // run simulation
    int num_eval_xcoords_per_element = 5;
    int num_eval_tcoords_per_element = 3;
	test.run(num_eval_xcoords_per_element,num_eval_tcoords_per_element,
             ufile,qfile,xfile,fufile,fqfile,fxfile);

    // print final message
	cout<<"\nThe Simulation has completed!\n";
	return 0;
}

