

#ifndef CELL_H
#define CELL_H

#include "vector.h"
#include "geometry.h"
#include "shapefunctions.h"
#include "calculus.h"
#include "matrix.h"


enum Field{Conc,Flux};

double initialCondition(double x );

class Cell
{
	public:
		//Constructors
		Cell():row(0),col(0),p(0),Dt(0),Dx(0){}
		Cell(int r, int c, int P,double Dt, double Dx);

		//Function to change the cell characteristics
		void change(int r, int c, int P, double dt, double dx);
		void change_pOrder(int P){ p = P; } 

		//Checks for where it is in the domain
		bool ifInitial(){ return row == 1; }
		bool ifFinal(int Nt){ return row == Nt; }
		bool ifLeftBoundary(){ return col == 1; }
		bool ifRightBoundary(int Nx){ return col == Nx; }


		//Functions for getting coordinate for its neighbors
		Point getRightNeighbor(int Nx);
		Point getLeftNeighbor(int Nx);
		Point getTopNeighbor(int Nt);
		Point getBottomNeighbor();

		//Check if a given coordinate in the global domain is within
		//the given element's domain
		bool withinCell( const Point & point) const;
		bool withinCell( double x, double y) const;

		//Evaluate the function at the coordinate coord
		double evaluateN(double zeta, double eta,int dx,int dy, int field,Vector const & Ncoeffs, Points const & Sf_arr);
		double evaluate(double x, double y, int field,Vector const & Ncoeffs, Points const & Sf_arr);
		double evaluate(double x, double y, int dx, int dy, int field,Vector const & Ncoeffs, Points const & Sf_arr);
		Vector evaluate(Vector x, double y, int field,Vector const & Ncoeffs, Points const & Sf_arr);
		Vector evaluate(Vector x, double y, int dx, int dt, int field,Vector const & Ncoeffs, Points const & Sf_arr);
		double evaluateShapeFunc(double zeta,double eta,int dz, int de,
			                     int which, Vector const & Ncoeffs, Points const & sf_arr) const;
	
		//Transform (x,y) coordinates into natural ones, (zeta,eta)
		double x2natural(double x) const;
		double y2natural(double y) const;
		double natural2x(double zeta) const;
		double natural2y(double eta) const;

		//Get the Jacobian value
		double getJacobian()const { return Jacobian; }

		//Get dx or dy
		double getDx()const { return Dx; }
		double getDy()const { return Dt; }
		double getP()const { return p; }

		//Get coefficients
		double getU(int i)const { return ucoeffs[i];}
		double getQ(int i)const { return qcoeffs[i];}
		double & getU(int i){ return ucoeffs[i];}
		double & getQ(int i){ return qcoeffs[i];}
	
		//Print location of cell in overall grid
		void print();

		//Print the info about the cell
		void printInfo();

	private:

		//Row and column location of this cell
		int row;
		int col;

		//Width of cell
		double Dx;

		//Height of cell
		double Dt;

		//Polynomial interpolation order
		int p;

		//Jacobian value for integration
		double Jacobian;

		//Concentration solution values
		Vector ucoeffs;

		//Flux solution values
		Vector qcoeffs;

		//Get coordinates of vertices
		Point node1() const;   //Bottom left point
		Point node2() const;   //Bottom right point
		Point node3() const;   //Top right point
		Point node4() const;   //Top left point

		//Get left x value
		double leftx() const { return (col-1)*Dx; }
		//Get right x value
		double rightx() const { return col*Dx; }
		//Get top y value
		double topy() const { return row*Dt; }
		//Get bottom y value
		double bottomy() const { return (row-1)*Dt; }
};


#endif