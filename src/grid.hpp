//
//  grid.hpp
//  TAM574_STDG
//
//  Created by Christian J Howard on 9/11/16.
//
//  The MIT License (MIT)
//    Copyright © 2016 Christian Howard. All rights reserved.
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in all
//  copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//  SOFTWARE.
//
//

#ifndef grid_hpp
#define grid_hpp

#include "vector.h"
#include "geometry.h"
#include "shapefunctions.h"
#include "calculus.h"
#include "matrix.h"
#include "cell.h"

class Grid{
    
public:
    Grid():grid(0),rows(0),cols(0){}
    Grid(int rows, int cols, int poly_order, double width, double height, const Vector & v);
    void change(int rows, int cols, int poly_order, double width, double height, const Vector & v);
    void assignConstants(Vector const & v){ constants = v; }
    ~Grid();
    Cell & operator()(int row, int col); //Using (1,1) as min (row,col)
    Grid & operator=( const Grid grid);
    Vector ncoeffs;//Normalized coefficients for weight functions
    Points sf_arr; //Shape function array. Will show [(0,0),(1,0),(0,1),...] which is the
    //polynomial orders for zeta and eta for each interpolation term
    Vector weights;
    Vector abscissa;
    
    double volumeIntegral(int r1,int c1,
                          int dz1,int de1,int dz2,int de2,
                          int which1,int which2) const;
    double lineIntegral(int r1,int c1,int r2,int c2,
                        int dx1,int dy1,int dx2,int dy2,
                        int which1,int which2, int edge,int field) const;
    double u_star(double param, int r, int c,int edge);
    double q_star(double param, int r,int c,int edge);
    double R_u(int r,int c);
    double Denom_u(int r,int c);
    double R_q(int r,int c);
    double Denom_q(int r,int c);
    
    double convergence_u();
    double convergence_q();
    
    //Local equations for stiffness matrix
    double am1(int currR,int currC,int which1,int which2) const;
    double a0(int currR,int currC,int which1,int which2) const;
    double ap1(int currR,int currC,int which1,int which2) const;
    
    double bm1(int currR,int currC,int which1,int which2) const;
    double b0(int currR,int currC,int which1,int which2) const;
    double bp1(int currR,int currC,int which1,int which2) const;
    
    double cm1(int currR,int currC,int which1,int which2) const;
    double c0(int currR,int currC,int which1,int which2) const;
    double cp1(int currR,int currC,int which1,int which2) const;
    
    double dm1(int currR,int currC,int which1,int which2) const;
    double d0(int currR,int currC,int which1,int which2) const;
    double dp1(int currR,int currC,int which1,int which2) const;
    
    double F(int currR,int currC, int which) const; //Load vector portion from weight function v
    double G(int currR,int currC, int which) const; //Load vector portion from weight function w
    
    void GenerateStiffnessMatrix(Matrix & K);
    void GenerateLoadVector(int row);
    void SolveConstants(int row);
    void getInverse();
    double (*initU)(double x);
    double (*initQ)(double x);
    
    bool empty(){ return (rows <1 || cols < 1);}
    int Rows(){ return rows; }
    int Cols(){ return cols; }
    
    void writeSolution2File(int numX,int numT,string fileu, string fileq,string filex);
    void writeEndSolution2File(int numx,int numT,string fileu,string fileq,string filex);
    //numX is number of evaluations per cell in the spacial dimension
    //numT is number of evaluations per cell in the time dimension
    
private:
    Cell** grid;
    Vector constants;
    Cell null;
    
    int rows;
    int cols;
    void clear();
    int e2gR(int weight_elem,int weight_field,int weight_func);
    int e2gC(int soln_elem,int soln_field,int soln_func);
    
    Matrix Inverse;
    
    Vector Load;
    
};

#endif /* grid_hpp */
