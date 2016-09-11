
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include "vector.h"
using namespace std;

class Matrix
{

public:

	Matrix():rows(0),cols(0),arr(NULL){};
	Matrix(int row, int col);
	Matrix(int row, int col,double all_value);
	Matrix(const Matrix & m);
	~Matrix();
	void resize(int r,int c);
	Matrix operator+(const Matrix & m) const;
	Matrix operator-(const Matrix & m) const;
	Matrix & operator=(const Matrix & m);
	Matrix operator*(const Matrix & m) const;
	Vector operator*(const Vector & v) const;
	double operator()(int r, int c) const;
	double & operator()(int r, int c);
	Matrix gaussSolve(Matrix b);
	Matrix gaussSolve(Vector b);
	Vector gaussSolveV(Vector b);
	void diag(int dim,double value);
	void invert(Matrix & m);
	void print();
	void write2file(string filename);

	int Rows() const{ return rows;}
	int Cols() const{ return cols;}

	void multrow(int row,double value);
	void swaprow(int row1,int row2);
	void addrow(int row_toadd,double bycoeff,int row_addto);

private:

	double** arr;
	int rows;
	int cols;
	void clear();
	void copy(const Matrix & m);
};


class Lower{
public:
	Lower():arr(0),dim(0){}
	Lower(int d);
	Lower(int d,double value);
	~Lower();
	void resize(int dim);
	void toDiagonal(double val);
	Vector solve(Vector const & b);
	double Rows(){ return dim;}
	double Cols(){ return dim;}
	double determinant();
	double operator()(int r, int c) const;
	double & operator()(int r, int c);
	void print();

private:
	double**arr;
	int dim;
	void clear();
	void copy( const Lower * l);

};


class Upper{
public:
	Upper():arr(0),dim(0){}
	Upper(int d);
	Upper(int d,double value);
	~Upper();
	void resize(int dim);
	void toDiagonal(double val);
	Vector solve(Vector const & b);
	double Rows(){ return dim;}
	double Cols(){ return dim;}
	double determinant();
	double operator()(int r, int c) const;
	double & operator()(int r, int c);
	void print();

private:
	double**arr;
	int dim;
	void clear();
	void copy( const Upper * l);

};


void LU_Decomposition(Lower & L, Upper & U, const Matrix & A);
#endif