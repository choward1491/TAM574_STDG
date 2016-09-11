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
	double T = 1;
	double L = 2*PI;
	double a = 1;
	consts[0] = 2*PI;	//a
	consts[1] = consts[0];	//c
	consts[2] = 1.0/(L*L);	//tau
	consts[3] = 1;	//kappa
	int Nx = 3;
	int Nt = 40;
	int m = 1;

	Simulation test(Nt,Nx,2,L,T,&FunctionU,&FunctionQ,consts);

	test.run(5,3,"u_p3_c1.txt","q_p3_c1.txt","x_p3_c1.txt");

	cout<<"\nThe Simulation has completed!\n";
	return 0;
}

