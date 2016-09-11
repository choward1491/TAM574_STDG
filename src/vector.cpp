
#include "vector.h"


Vector::Vector(int size)
{
	len = size;
	arr = new double[size]();
}

void Vector::change(int size)
{
	clear();
	len = size;
	arr = new double[size]();
}

Vector::Vector( const Vector & v)
{
	copy(v);
}

Vector::Vector(double* & arr, int size)
{
	len = size;
	arr = new double[size];

	for(int i = 0; i < len; i++)
		arr[i] = arr[i];
}

void Vector::print()
{
	if(arr != NULL)
	{
		cout<<"<";
		for(int i = 0; i < len; i++)
			cout<<" "<<arr[i];
		cout<<" >\n";
	}
}

void Vector::clear()
{
	if( arr != NULL )
	{
		delete [] arr;
		arr = NULL;
	}
	len = 0;
}

void Vector::copy(const Vector & v)
{
	arr = new double[v.length()]();
	len = v.length();

	for(int i = 0; i < len; i++)
		arr[i] = v[i];
}

Vector::~Vector()
{
	if( arr != NULL)
		clear();
}

Vector::Vector(double start, double end, int numpts)
{
	len = numpts;
	arr = new double[numpts];

	double dx = (end-start)/(double)(numpts-1);

	for( int i = 0; i<numpts; i++)
		arr[i] = start+i*dx;
}


Vector Vector::operator+(const Vector& v2) const
{
	if( v2.length() == len )
	{
		Vector out(len);

		for(int i = 0; i < len; i++ ){
			out[i] = arr[i] + v2[i];
		}

		return out;
	}

	cout<<"arr sizes don't match, returning empty vector\n";
	return Vector();
}
Vector Vector::operator-(const Vector & v2) const
{
	if( v2.length() == len )
	{
		Vector out(len);

		for(int i = 0; i < len; i++ ){
			out[i] = arr[i] - v2[i];
		}

		return out;
	}

	cout<<"arr sizes don't match, returning empty vector\n";
	return Vector();
}

Vector & Vector::operator=(const Vector & v)
{
	if( &v != this)
	{
		clear();
		copy(v);
	}

	return *this;
}

void Vector::write2file(string filename)
{
	if( arr != 0 )
	{

		ofstream Fu;
		Fu.open(filename,ios::app);

		//Do file for u first
		for( int i = 0; i < len; i++)
		{
					Fu<<arr[i]<<" ";
		}Fu<<"\n";
		Fu.close();

	}
}