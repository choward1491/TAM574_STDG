
#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

class Vector
{

public:
	Vector():len(0),arr(NULL){}
	Vector( const Vector & v);
	Vector(int size);
	Vector(double*& Array, int size);
	Vector(double start, double end, int numpts);
	~Vector();
	int length() const { return len;}
	double operator[](int index)const { return arr[index];}
	double & operator[](int index){ return arr[index];}
	void change(int size);
	Vector operator+(const Vector & v2) const;
	Vector operator-(const Vector & v2) const;
	Vector & operator=(const Vector & v);
	void write2file(string filename);
	void print();
	
private:
	int len;
	double* arr;
	void clear();
	void copy(const Vector & v);
};

#endif