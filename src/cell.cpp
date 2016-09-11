

#include "cell.h"




Cell::Cell(int r, int c, int P,double dt, double dx)
{
	row = r; col = c; p = P; Dt = dt; Dx = dx;
	int tmp = simplesum(1,p+1,1);
	ucoeffs.change(tmp);
	qcoeffs.change(tmp);

	Jacobian = dx*dt/4.0;
}

void Cell::change(int r, int c, int P, double dt, double dx)
{
	row = r; col = c; p = P; Dt = dt; Dx = dx;
	int tmp = simplesum(1,p+1,1);
	ucoeffs.change(tmp);
	qcoeffs.change(tmp);
	Jacobian = dx*dt/4.0;
}

Point Cell::getRightNeighbor(int Nx)
{
	Point out;

	if( ifRightBoundary(Nx) )
	{
		out.change(row,1);
		return out;
	}
	
	out.change(row,col+1);

	return out;
}
Point Cell::getLeftNeighbor(int Nx)
{
	Point out;

	if( ifLeftBoundary() )
	{
		out.change(row,Nx);
		return out;
	}
	
	out.change(row,col-1);
	return out;
}
Point Cell::getTopNeighbor(int Nt)
{
	Point out;

	if( ifFinal(Nt) )
	{
		out.change(-1,-1);
		return out;
	}
	
	out.change(row+1,col);
	return out;
}
Point Cell::getBottomNeighbor()
{
	Point out;

	if( ifInitial() )
	{
		out.change(-1,-1);
		return out;
	}
	
	out.change(row-1,col);
	return out;
}

bool Cell::withinCell( const Point & point) const
{
	double x = point.x, y = point.y;
	double xL = (col-1)*Dx, xR = col*Dx;
	double yB = Dt*(row-1), yT = row*Dt;

	if( x <= xR && x >= xL && y <= yT && y >= yB)
		return true;
	return false;
}

bool Cell::withinCell( double x, double y) const
{
	double xL = (col-1)*Dx, xR = col*Dx;
	double yB = Dt*(row-1), yT = row*Dt;

	if( x <= xR && x >= xL && y <= yT && y >= yB)
		return true;
	return false;
}

Point Cell::node1() const   //Bottom left point
{
	Point out;
	out.x = (col-1)*Dx;
	out.y = (row-1)*Dt;
	return out;
}
Point Cell::node2() const   //Bottom right point
{
	Point out;
	out.x = (col)*Dx;
	out.y = (row-1)*Dt;
	return out;
}
Point Cell::node3() const   //Top right point
{
	Point out;
	out.x = (col)*Dx;
	out.y = (row)*Dt;
	return out;
}
Point Cell::node4() const   //Top left point
{
	Point out;
	out.x = (col-1)*Dx;
	out.y = (row)*Dt;
	return out;
}

double Cell::x2natural(double x) const
{
	double xL = leftx(), xR = rightx();
	return (2*x - (xL+xR))/Dx;
}
double Cell::y2natural(double y) const
{
	double yB = bottomy(), yT = topy();
	return (2*y - (yB+yT))/Dt;
}

double Cell::natural2x(double zeta) const
{
	double xL = leftx(), xR = rightx();
	return (zeta*Dx + (xL+xR) )/2.0;
}
double Cell::natural2y(double eta) const
{
	double yB = bottomy(), yT = topy();
	return (eta*Dt + (yT+yB) )/2.0;
}

double Cell::evaluateN(double zeta, double eta, int dx,int dy,int field,Vector const & Ncoeffs, Points const & Sf_arr)
{
	double output = 0;

	Vector* ptr = 0;
	Point* pt = 0;

	if( field == 0 )
		ptr = &ucoeffs;
	else
		ptr = &qcoeffs;

	for(int i = 0; i < Ncoeffs.length(); i++)
	{
		output += (*ptr)[i]*evaluateShapeFunc(zeta,eta,dx,dy,i,Ncoeffs,Sf_arr);
	}
	return output;
}

double Cell::evaluate(double x, double y, int field,Vector const & Ncoeffs, Points const & Sf_arr)
{
	if( withinCell(x,y) )
	{
		Vector* tmp = 0;

		if(field == 0)
			tmp = &ucoeffs;
		else
			tmp = &qcoeffs;

		int len = Ncoeffs.length();
		double zeta = x2natural(x);
		double eta = y2natural(y);
		double out = 0;

		for(int i = 0; i < len; i++)
		{
			out += (*tmp)[i]*evaluateShapeFunc(zeta,eta,0,0,i,Ncoeffs,Sf_arr);
		}

		return out;
	}

	return 0;
}

double Cell::evaluate(double x, double y, int dx, int dy, int field,Vector const & Ncoeffs, Points const & Sf_arr)
{
	if( withinCell(x,y) )
	{
		Vector* tmp = 0;

		if(field == 0)
			tmp = &ucoeffs;
		else
			tmp = &qcoeffs;

		int len = Ncoeffs.length();
		double zeta = x2natural(x);
		double eta = y2natural(y);
		double out = 0;

		for(int i = 0; i < len; i++)
		{
			out += (*tmp)[i]*evaluateShapeFunc(zeta,eta,dx,dy,i,Ncoeffs,Sf_arr);
		}

		return out;
	}

	return 0;
}

Vector Cell::evaluate(Vector x, double y, int field,Vector const & Ncoeffs, Points const & Sf_arr)
{
	Vector out(x.length());

	for(int i = 0; i < x.length(); i++)
	{
		out[i] = evaluate(x[i], y, field, Ncoeffs, Sf_arr);
	}

	return out;
}

Vector Cell::evaluate(Vector x, double y, int dx, int dy, int field,Vector const & Ncoeffs, Points const & Sf_arr)
{
	Vector out(x.length());
	double tmp = 0;

	for(int i = 0; i < x.length(); i++)
	{
		out[i] = evaluate(x[i],y,dx,dy,field,Ncoeffs,Sf_arr);
	}

	return out;
}

double Cell::evaluateShapeFunc(double zeta,double eta,
	                           int dx, int dy,
							   int which, //Number for the shape function in vector
							   Vector const & Ncoeffs, 
							   Points const & Sf_arr) const //Finds the value of shape function, including
														    //derivatives of the real coordinates
{
	Point pt = Sf_arr[which];
	double coeff = Ncoeffs[which];

	if(dx > pt.x || dy > pt.y ) //If the desired derivative is more than the order of polynomial
		return 0;
	else{                       //Otherwise, obtain the nonzero value
		double zord = pt.x-dx, eord = pt.y - dy;
		double c1 = 1,c2 = 1;

		for( int i = pt.x; i > zord; i--)
			c1*=i;
		for( int i = pt.y; i > eord; i--)
			c2*=i;

		return c1*c2*coeff*pow(zeta,zord)*pow(eta,eord)*pow(2.0/Dx,dx)*pow(2.0/Dt,dy);
	}
}

void Cell::print()
{
	cout<<"("<<row<<", "<<col<<") ";
}

void Cell::printInfo()
{
	cout<<"This cell is located at "; print(); cout<<"\n";
	cout<<"This cell has complete interpolant polynomials \n  with order "<<p<<" interpolation.\n";
	cout<<"This cell has a width of "<<Dx<<" and a height of "<<Dt<<".\n";
}
