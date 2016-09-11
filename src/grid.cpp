//
//  grid.cpp
//  TAM574_STDG
//
//  Created by Christian J Howard on 9/11/16.
//  Copyright © 2016 Christian Howard. All rights reserved.
//

#include "grid.hpp"

void Grid::clear()
{
    if( grid != 0 )
    {
        for( int i = 0; i < rows; i++)
            delete [] grid[i];
        delete [] grid;
    }
    grid = 0;
}

Grid::~Grid()
{
    clear();
}

Grid::Grid(int r, int c, int poly_order, double width, double height, const Vector & v)
{
    //Get the rows and columns of the grid
    rows = r; cols = c;
    
    //Obtain constants for simulation
    assignConstants(v);
    
    //Obtain dx and dy for each element
    double dx = width/c;
    double dy = height/r;
    
    //Number of terms for interpolants
    int tmp = simplesum(1,poly_order+1,1);
    
    //Get normalized coefficients for shape functions
    ncoeffs = getNormedShapeCoeffs(poly_order,dx,dy);
    
    //Get the order of each natural variable for each coefficient
    sf_arr = getShapeFuncPolyOrders(poly_order);
    
    //Get the weights and abscissa
    Vector quad = getQuadrature(poly_order);
    
    //Separate the quadrature schemes into weights and abscissa vector
    int len = quad.length();
    weights.change(len/2);
    abscissa.change(len/2);
    
    for(int i = 0; i < len/2; i++)
    {
        abscissa[i] = quad[i];
        weights[i] = quad[len/2+i];
    }
    
    //Create the grid
    grid = new Cell*[rows];
    
    for(int i = 0; i < rows; i++ )
        grid[i] = new Cell[cols];
    
    for( int i = 0; i < rows; i++)
    {
        for( int j = 0; j < cols; j++ )
            grid[i][j].change(i+1,j+1,poly_order,dy,dx);
    }
    
    //Get the Inverse of stiffness matrix
    getInverse();
    
    //End
}

void Grid::change(int r, int c, int poly_order, double width, double height, const Vector & v)
{
    clear();
    
    //Get the rows and columns of the grid
    rows = r; cols = c;
    
    //Obtain constants for simulation
    assignConstants(v);
    
    //Obtain dx and dy for each element
    double dx = width/c;
    double dy = height/r;
    
    //Number of terms for interpolants
    int tmp = simplesum(1,poly_order+1,1);
    
    //Get normalized coefficients for shape functions
    ncoeffs = getNormedShapeCoeffs(poly_order,dx,dy);
    
    //Get the order of each natural variable for each coefficient
    sf_arr = getShapeFuncPolyOrders(poly_order);
    
    //Get the weights and abscissa
    Vector quad = getQuadrature(poly_order);
    
    //Separate the quadrature schemes into weights and abscissa vector
    int len = quad.length();
    weights.change(len/2);
    abscissa.change(len/2);
    
    for(int i = 0; i < len/2; i++)
    {
        abscissa[i] = quad[i];
        weights[i] = quad[len/2+i];
    }
    
    //Create the grid
    grid = new Cell*[rows];
    
    for(int i = 0; i < rows; i++ )
        grid[i] = new Cell[cols];
    
    for( int i = 0; i < rows; i++)
    {
        for( int j = 0; j < cols; j++ )
            grid[i][j].change(i+1,j+1,poly_order,dy,dx);
    }
    
    //Get the Inverse of stiffness matrix
    getInverse();
    
    //End
}


Cell & Grid::operator()(int row, int col)
{
    
    if( row > 0 && row <= rows)
    {
        if( col > 0 && col <= cols)
            return grid[row-1][col-1];
    }
    
    cout<<"Tried accessing cell ("<<row<<", "<<col<<") when the dims\n";
    cout<<" are [0,"<<rows<<"] x [0,"<<cols<<"]. Try again\n\n";
    null = Cell();
    return null;
}

double initialCondition(double x )
{
    return 0;
}

double Grid::volumeIntegral(int r1,int c1,
                            int dx1,int dy1,int dx2,int dy2,
                            int which1,int which2) const
{
    double output = 0;
    
    if( grid != 0 )
    {
        
        Cell *pt1 = &grid[r1][c1];
        
        int len = weights.length();
        
        for( int i = 0; i < len; i++)
        {
            for( int j = 0; j < len; j++ )
            {
                
                output += weights[i]*weights[j]*
                pt1->evaluateShapeFunc(abscissa[i],abscissa[j],dx1,dy1,which1, ncoeffs, sf_arr)*
                pt1->evaluateShapeFunc(abscissa[i],abscissa[j],dx2,dy2,which2, ncoeffs, sf_arr);
            }
        }
        
        output*=pt1->getJacobian();
        
    }
    
    return output;
}

double Grid::lineIntegral(int r1,int c1,int r2,int c2,
                          int dx1,int dy1,int dx2,int dy2,
                          int which1,int which2,int edge,int field) const
{
    double output = 0;
    
    //initial jacobian value
    double jacobian = .5;
    
    if( grid != 0 )
    {
        Cell *pt1 = &grid[r1][c1], *pt2 = &grid[r2][c2];
        int len = weights.length();
        double const1 = 0,const2 = 0;
        
        switch(edge)
        {
            case 1:
                jacobian *= pt1->getDx();
                const1 = -1;
                if( r1 == 0 && r2 == -1)
                {	double X = 0;
                    double W, sf,IC;
                    if( field == 0 )
                    {
                        for( int i = 0; i < len; i++)
                        {
                            X = pt1->natural2x(abscissa[i]);
                            W = weights[i];
                            sf = pt1->evaluateShapeFunc(abscissa[i],const1,dx1,dy1,which1, ncoeffs, sf_arr);
                            IC = initU(X);
                            output += W*IC*sf;
                        }
                    }else{
                        for( int i = 0; i < len; i++)
                        {
                            X = pt1->natural2x(abscissa[i]);
                            output += weights[i]*
                            pt1->evaluateShapeFunc(abscissa[i],const1,dx1,dy1,which1, ncoeffs, sf_arr)*initQ(X);
                        }
                    }
                    
                }else{
                    if( r2 != r1 ){
                        const2 = -const1;
                    }else
                        const2 = const1;
                    
                    for( int i = 0; i < len; i++)
                    {
                        if( r1 == r2)
                        {
                            output += weights[i]*
                            pt1->evaluateShapeFunc(abscissa[i],const1,dx1,dy1,which1, ncoeffs, sf_arr)*
                            pt2->evaluateShapeFunc(abscissa[i],const2,dx2,dy2,which2, ncoeffs, sf_arr);
                        }else{
                            output += weights[i]*
                            pt1->evaluateShapeFunc(abscissa[i],const1,dx1,dy1,which1, ncoeffs, sf_arr)*
                            pt2->evaluateN(abscissa[i], const2,0,0, field,ncoeffs, sf_arr);
                        }
                    }
                }
                break;
            case 2:
                jacobian *= pt1->getDy();
                const1 = 1;
                if( c2 != c1 ){
                    const2 = -const1;
                }else
                    const2 = const1;
                
                for( int i = 0; i < len; i++)
                {
                    output += weights[i]*
                    pt1->evaluateShapeFunc(const1,abscissa[i],dx1,dy1,which1, ncoeffs, sf_arr)*
                    pt2->evaluateShapeFunc(const2,abscissa[i],dx2,dy2,which2, ncoeffs, sf_arr);
                }
                
                break;
                
            case 3:
                jacobian *= -pt1->getDx();
                const1 = 1;
                
                if( r2 != r1 ){
                    const2 = -const1;
                }else
                    const2 = const1;
                
                for( int i = 0; i < len; i++)
                {
                    output += weights[i]*
                    pt1->evaluateShapeFunc(abscissa[i],const1,dx1,dy1,which1, ncoeffs, sf_arr)*
                    pt2->evaluateShapeFunc(abscissa[i],const2,dx2,dy2,which2, ncoeffs, sf_arr);
                }
                break;
            case 4:
                jacobian *= -pt1->getDy();
                const1 = -1;
                if( c2 != c1 ){
                    const2 = -const1;
                }else
                    const2 = const1;
                
                for( int i = 0; i < len; i++)
                {
                    output += weights[i]*
                    pt1->evaluateShapeFunc(const1,abscissa[i],dx1,dy1,which1, ncoeffs, sf_arr)*
                    pt2->evaluateShapeFunc(const2,abscissa[i],dx2,dy2,which2, ncoeffs, sf_arr);
                }
                
                break;
            default:
                cout<< "No edge value within the range of 1-4 was used, so returning 0\n";
                break;
        }
        
        output *= jacobian;
    }
    
    return output;
}

/*
 *
 *
 *
 *
 *     Obtain the local stiffness stuff!
 *
 *
 *
 *     Note that for this situation, the following is true:
 *
 *	  constants = { a, c, tau, k }
 *
 *
 *
 *
 */


/*
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *     This is the a part of the equations
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 */

double Grid::am1(int currR,int currC,int which1,int which2) const
{
    double tau = constants[2];
    double k = constants[3];
    double a = constants[0];
    double c = constants[1];
    int tmp = 0;
    double C = (a+c)/2.0;
    
    if( currC == 0 )
        tmp = cols-1;
    else
        tmp = currC-1;
    
    //Boundary integral components
    return C* lineIntegral(currR,currC,currR,tmp,0,0,0,0,which1,which2,4,0);
}

double Grid::a0(int currR,int currC,int which1,int which2) const
{
    double output = 0;
    
    double tau = constants[2];
    double k = constants[3];
    double a = constants[0];
    double c = constants[1];
    
    double C = 0;
    
    //Volume integral components
    output += -volumeIntegral(currR,currC,0,0,0,1,which1,which2);
    output += -a*volumeIntegral(currR,currC,0,0,1,0,which1,which2);
    
    //Boundary integral components
    output += .5*(a+c)*lineIntegral(currR,currC,currR,currC,0,0,0,0,which1,which2, 2,0);
    output += .5*(a-c)*lineIntegral(currR,currC,currR,currC,0,0,0,0,which1,which2, 4,0);
    output += -lineIntegral(currR,currC,currR,currC,0,0,0,0,which1,which2, 3,0);
    
    //Return the output
    return output;
}

double Grid::ap1(int currR,int currC,int which1,int which2) const
{
    double tau = constants[2];
    double k = constants[3];
    double a = constants[0];
    double c = constants[1];
    int tmp = 0;
    double C = (a-c)/2.0;
    
    if( currC == cols-1 )
        tmp = 0;
    else
        tmp = currC+1;
    
    //Boundary integral components
    return C* lineIntegral(currR,currC,currR,tmp,0,0,0,0,which1,which2,2,0);
}


/*
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *     This is the b part of the equations
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 */
double Grid::bm1(int currR,int currC,int which1,int which2) const
{
    double tau = constants[2];
    double k = constants[3];
    double a = constants[0];
    double c = constants[1];
    int tmp = 0;
    double C = (a/c+1)/2.0;
    
    if( currC == 0 )
        tmp = cols-1;
    else
        tmp = currC-1;
    
    //Boundary integral components
    return C* lineIntegral(currR,currC,currR,tmp,0,0,0,0,which1,which2,4,0);
}

double Grid::b0(int currR,int currC,int which1,int which2) const
{
    double output = 0;
    
    double tau = constants[2];
    double k = constants[3];
    double a = constants[0];
    double c = constants[1];
    
    double C = 0;
    
    //Volume integral components
    output += -volumeIntegral(currR,currC,0,0,1,0,which1,which2);
    
    //Boundary integral components
    output += .5*(1+a/c)*lineIntegral(currR,currC,currR,currC,0,0,0,0,which1,which2, 2,0);
    output += .5*(1-a/c)*lineIntegral(currR,currC,currR,currC,0,0,0,0,which1,which2, 4,0);
    
    //Return output
    return output;
}

double Grid::bp1(int currR,int currC,int which1,int which2) const
{
    double tau = constants[2];
    double k = constants[3];
    double a = constants[0];
    double c = constants[1];
    int tmp = 0;
    double C = (1-a/c)/2.0;
    
    if( currC == cols-1 )
        tmp = 0;
    else
        tmp = currC+1;
    
    //Boundary integral components
    return C* lineIntegral(currR,currC,currR,tmp,0,0,0,0,which1,which2,2,0);
}


/*
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *     This is the c part of the equations
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 */
double Grid::cm1(int currR,int currC,int which1,int which2) const
{
    double tau = constants[2];
    double k = constants[3];
    double a = constants[0];
    double c = constants[1];
    int tmp = 0;
    double C = (k+a*tau*c)/2.0;
    
    if( currC == 0 )
        tmp = cols-1;
    else
        tmp = currC-1;
    
    //Boundary integral components
    return C* lineIntegral(currR,currC,currR,tmp,0,0,0,0,which1,which2,4,0);
}

double Grid::c0(int currR,int currC,int which1,int which2) const
{
    double output = 0;
    
    double tau = constants[2];
    double k = constants[3];
    double a = constants[0];
    double c = constants[1];
    
    double C = 0;
    
    //The volume integral components
    output += -k*volumeIntegral(currR,currC,0,0,1,0,which1,which2);
    
    //The boundary integral components
    output += .5*(k+a*tau*c)*lineIntegral(currR,currC,currR,currC,0,0,0,0,which1,which2, 2,0);
    output += .5*(k-a*tau*c)*lineIntegral(currR,currC,currR,currC,0,0,0,0,which1,which2, 4,0);
    
    //Return the value
    return output;
}

double Grid::cp1(int currR,int currC,int which1,int which2) const
{
    double tau = constants[2];
    double k = constants[3];
    double a = constants[0];
    double c = constants[1];
    int tmp = 0;
    double C = (k-a*tau*c)/2.0;
    
    if( currC == cols-1 )
        tmp = 0;
    else
        tmp = currC+1;
    
    return C* lineIntegral(currR,currC,currR,tmp,0,0,0,0,which1,which2,2,0);
}


/*
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *     This is the d part of the equations
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 */
double Grid::dm1(int currR,int currC,int which1,int which2) const
{
    double tau = constants[2];
    double k = constants[3];
    double a = constants[0];
    double c = constants[1];
    int tmp = 0;
    double C = (k/c + a*tau)/2.0;
    
    if( currC == 0 )
        tmp = cols-1;
    else
        tmp = currC-1;
    
    return C* lineIntegral(currR,currC,currR,tmp,0,0,0,0,which1,which2,4,0);
}

double Grid::d0(int currR,int currC,int which1,int which2) const
{
    double output = 0;
    
    double tau = constants[2];
    double k = constants[3];
    double a = constants[0];
    double c = constants[1];
    
    double C = 0;
    
    //Volume integral terms
    output += volumeIntegral(currR,currC,0,0,0,0,which1,which2);
    output += -tau*volumeIntegral(currR,currC,0,0,0,1,which1,which2);
    output += -a*tau*volumeIntegral(currR,currC,0,0,1,0,which1,which2);
    
    //Boundary integral terms
    output += .5*(k/c + a*tau)*lineIntegral(currR,currC,currR,currC,0,0,0,0,which1,which2, 2,0);
    output += .5*(a*tau-k/c)*lineIntegral(currR,currC,currR,currC,0,0,0,0,which1,which2, 4,0);
    output += -tau*lineIntegral(currR,currC,currR,currC,0,0,0,0,which1,which2, 3,0);
    
    //Return the output
    return output;
}

double Grid::dp1(int currR,int currC,int which1,int which2) const
{
    double tau = constants[2];
    double k = constants[3];
    double a = constants[0];
    double c = constants[1];
    int tmp = 0;
    double C = (a*tau-k/c)/2.0;
    
    if( currC == cols-1 )
        tmp = 0;
    else
        tmp = currC+1;
    
    return C* lineIntegral(currR,currC,currR,tmp,0,0,0,0,which1,which2,2,0);
}


/*
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *     This is the loading part of the equations
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 */
double Grid::F(int currR,int currC, int which) const //Load vector portion from weight function v
{
    double output = 0;
    
    output = lineIntegral(currR,currC,currR-1,currC,0,0,0,0,which,which,1,0);
    
    return output;
}
double Grid::G(int currR,int currC, int which) const //Load vector portion from weight function w
{
    double output = 0;
    
    double tau = constants[2];
    double k = constants[3];
    double a = constants[0];
    double c = constants[1];
    
    output = tau*lineIntegral(currR,currC,currR-1,currC,0,0,0,0,which,which,1,1);
    return output;
}


/*
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *     The formal compilation of the equations into the
 *                  global stiffness matrix
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 */

int Grid::e2gR(int weight_elem,int weight_field,int weight_func)
{
    int out;
    int num_p = ncoeffs.length();
    int num_elems = cols;
    
    return (weight_elem+weight_field*num_elems)*num_p+weight_func;
}
int Grid::e2gC(int soln_elem,int soln_field,int soln_func)
{
    int out;
    int num_p = ncoeffs.length();
    int num_elems = cols;
    
    return (soln_elem+soln_field*num_elems)*num_p+soln_func;
}

void Grid::GenerateStiffnessMatrix( Matrix & K)
{
    int num_b = ncoeffs.length();                 //number of basis functions per field per element
    int should_be = 2*num_b*cols;                    //how large the matrix K should be
    int minus,plus;
    
    if( K.Cols() != should_be || K.Rows() != should_be )
        K.resize(should_be,should_be);
    
    //Loop through elements
    for(int e = 0; e < cols; e++)
    {
        if( e == 0 )
            minus = cols-1;
        else
            minus = e-1;
        
        if(e == cols-1)
            plus = 0;
        else
            plus = e+1;
        
        //Loop through basis functions
        for(int n = 0; n < num_b; n++) //weight functions
        {
            for(int m = 0; m < num_b; m++) //function basis functions
            {
                
                K(e2gR(e,0,n),e2gC(e,0,m)) += a0(0,e,m,n);
                K(e2gR(e,0,n),e2gC(minus,0,m)) += am1(0,e,n,m);
                K(e2gR(e,0,n),e2gC(plus,0,m)) += ap1(0,e,n,m);
                //cout<<"K(r,c) = "<<ap1(0,e,n,m)<<",G(r,c) = "<<"("<<e2gR(e,0,n)<<","<<e2gC(plus,0,m)<<"), (r,c) = ("<<n<<", "<<m<<"), element "<<e<<endl;
                
                K(e2gR(e,0,n),e2gC(e,1,m)) += b0(0,e,m,n);
                K(e2gR(e,0,n),e2gC(minus,1,m)) += bm1(0,e,n,m);
                K(e2gR(e,0,n),e2gC(plus,1,m)) += bp1(0,e,n,m);
                
                K(e2gR(e,1,n),e2gC(e,0,m)) += c0(0,e,m,n);
                K(e2gR(e,1,n),e2gC(minus,0,m)) += cm1(0,e,n,m);
                K(e2gR(e,1,n),e2gC(plus,0,m)) += cp1(0,e,n,m);
                
                K(e2gR(e,1,n),e2gC(e,1,m)) += d0(0,e,m,n);
                K(e2gR(e,1,n),e2gC(minus,1,m)) += dm1(0,e,n,m);
                K(e2gR(e,1,n),e2gC(plus,1,m)) += dp1(0,e,n,m);
                
            }//End of m loop
        }//End of n loop
    }//End of e loop
    
    //End of the stiffness matrix creation
    K.write2file("stiffness_matrix_p2Nx3.txt");
}

void Grid::getInverse()
{
    if( !empty() )
    {
        Matrix K;
        GenerateStiffnessMatrix(K);
        K.invert(Inverse);
        Inverse.write2file("theInverse_p2Nx3.txt");
    }
}

void Grid::GenerateLoadVector(int row)
{
    if( Load.length() != Inverse.Rows() )
        Load.change(Inverse.Rows());
    
    int num_p = ncoeffs.length();
    int num = cols*num_p;
    
    for(int e = 0; e < cols; e++)
    {
        for(int i = 0; i < num_p ; i++)
        {
            Load[e*num_p+i] = F(row,e,i);
            Load[i+e*num_p+num] = G(row,e,i);
        }
    }
    
    Load.write2file("theLoad_p2Nx3.txt");
}


void Grid::SolveConstants(int row)
{
    GenerateLoadVector(row);
    Vector y = Inverse*Load;
    y.write2file("solutioncoeffs_p2Nx3.txt");
    Cell* pt = 0;
    int num_p = ncoeffs.length();
    
    for(int e = 0; e < cols; e++)
    {
        pt = &grid[row][e];
        for(int n = 0; n < num_p; n++)
        {
            pt->getU(n) = y[e*num_p+n];
            pt->getQ(n) = y[(e+cols)*num_p+n];
        }
    }
    
    
}


void Grid::writeSolution2File(int numX,int numT,string fileu, string fileq,string filex)
{
    if( grid != 0 )
    {
        Vector x(-1,1,numX);
        Vector t(-1,1,numT);
        
        double dt = grid[0][0].getDy()/(numT-1);
        double dx = grid[0][0].getDx()/(numX-1);
        double currt = 0;
        
        ofstream Fu;
        Fu.open(fileu,ios::trunc);
        
        //Do file for u first
        for( int i = 0; i < rows; i++)
        {
            for(int k = (i>0); k < t.length(); k++)
            {
                Fu<<currt<<" ";
                for( int j = 0; j < cols; j++ )
                {
                    for( int l = (j>0); l < x.length(); l++ )
                    {
                        Fu<<grid[i][j].evaluateN(x[l],t[k],0,0,0,ncoeffs,sf_arr)<<" ";
                    }
                }
                Fu<<"\n";currt+=dt;
            }
        }
        Fu.close();
        
        Fu.open(fileq,ios::trunc);
        currt = 0;
        //Get file made for q next
        for( int i = 0; i < rows; i++)
        {
            for(int k = (i>0); k < t.length(); k++)
            {
                Fu<<currt<<" ";
                for( int j = 0; j < cols; j++ )
                {
                    for( int l = (j>0); l < x.length(); l++ )
                    {
                        Fu<<grid[i][j].evaluateN(x[l],t[k],0,0,1,ncoeffs,sf_arr)<<" ";
                    }
                }
                Fu<<"\n";currt+=dt;
            }
        }
        
        Fu.close();
        
        currt = 0;
        Fu.open(filex,ios::trunc);
        for( int j = 0; j < cols; j++ )
        {
            for( int l = (j>0); l < x.length(); l++ )
            {
                Fu<<currt<<" "; currt+=dx;
            }
        }
        Fu<<"\n";
        Fu.close();
    }
}

void Grid::writeEndSolution2File(int numX,int numT,string fileu,string fileq,string filex)
{
    if( grid != 0 )
    {
        Vector x(-1,1,numX);
        
        double dt = grid[0][0].getDy()/(numT-1);
        double dx = grid[0][0].getDx()/(numX-1);
        double currt = rows*dt;
        
        ofstream Fu;
        Fu.open(fileu,ios::trunc);
        
        //Do file for u first
        int i = rows-1;
        Fu<<currt<<" ";
        for( int j = 0; j < cols; j++ )
        {
            for( int l = (j>0); l < x.length(); l++ )
            {
                Fu<<grid[i][j].evaluateN(x[l],1,0,0,0,ncoeffs,sf_arr)<<" ";
            }
        }
        Fu<<"\n";
        
        Fu.close();
        
        Fu.open(fileq,ios::trunc);
        //Get file made for q next
        
        Fu<<currt<<" ";
        for( int j = 0; j < cols; j++ )
        {
            for( int l = (j>0); l < x.length(); l++ )
            {
                Fu<<grid[i][j].evaluateN(x[l],1,0,0,1,ncoeffs,sf_arr)<<" ";
            }
        }
        Fu<<"\n";
        
        Fu.close();
        
        currt = 0;
        Fu.open(filex,ios::trunc);
        for( int j = 0; j < cols; j++ )
        {
            for( int l = (j>0); l < x.length(); l++ )
            {
                Fu<<currt<<" "; currt+=dx;
            }
        }
        Fu<<"\n";
        Fu.close();
    }
}
double Grid::u_star(double param, int r, int c,int edge)
{
    Cell* cell = &grid[r][c];
    
    double um = 0;
    double up = 0;
    double qp = 0;
    double qm = 0;
    
    switch(edge)
    {
        case 1:
            if( r == 0 )
            {
                double x = cell->natural2x(param);
                return initU(x);
            }else{
                return grid[r-1][c].evaluateN(param,1,0,0,0,ncoeffs,sf_arr);
            }
            
            break;
        case 2:
            if( c == cols-1)
                c = 0;
            else
                c++;
            up = grid[r][c].evaluateN(-1,param,0,0,0,ncoeffs,sf_arr);
            um = cell->evaluateN(1,param,0,0,0,ncoeffs,sf_arr);
            qp = grid[r][c].evaluateN(-1,param,0,0,1,ncoeffs,sf_arr);
            qm = cell->evaluateN(1,param,0,0,1,ncoeffs,sf_arr);
            return (up+um)/2.0 + (qm-qp)/(2*constants[1]);
            
            break;
        case 3:
            return cell->evaluateN(param,1,0,0,0,ncoeffs,sf_arr);
            break;
        case 4:
            if( c == 0)
                c = cols-1;
            else
                c--;
            um = grid[r][c].evaluateN(1,param,0,0,0,ncoeffs,sf_arr);
            up = cell->evaluateN(-1,param,0,0,0,ncoeffs,sf_arr);
            qm = grid[r][c].evaluateN(1,param,0,0,1,ncoeffs,sf_arr);
            qp = cell->evaluateN(-1,param,0,0,1,ncoeffs,sf_arr);
            return (up+um)/2.0 + (qm-qp)/(2*constants[1]);
            
            break;
    }
    cell = 0;
    return 0;
}
double Grid::q_star(double param, int r,int c,int edge)
{
    Cell* cell = &grid[r][c];
    
    double um = 0;
    double up = 0;
    double qp = 0;
    double qm = 0;
    
    switch(edge)
    {
        case 1:
            if( r == 0 )
            {
                double x = cell->natural2x(param);
                return initQ(x);
            }else{
                return grid[r-1][c].evaluateN(param,1,0,0,1,ncoeffs,sf_arr);
            }
            
            break;
        case 2:
            if( c == cols-1)
                c = 0;
            else
                c++;
            up = grid[r][c].evaluateN(-1,param,0,0,0,ncoeffs,sf_arr);
            um = cell->evaluateN(1,param,0,0,0,ncoeffs,sf_arr);
            qp = grid[r][c].evaluateN(-1,param,0,0,1,ncoeffs,sf_arr);
            qm = cell->evaluateN(1,param,0,0,1,ncoeffs,sf_arr);
            return (qp+qm)/2.0 + constants[1]*(um-up)/(2.0);
            
            break;
        case 3:
            return cell->evaluateN(param,1,0,0,1,ncoeffs,sf_arr);
            break;
        case 4:
            if( c == 0)
                c = cols-1;
            else
                c--;
            um = grid[r][c].evaluateN(1,param,0,0,0,ncoeffs,sf_arr);
            up = cell->evaluateN(-1,param,0,0,0,ncoeffs,sf_arr);
            qm = grid[r][c].evaluateN(1,param,0,0,1,ncoeffs,sf_arr);
            qp = cell->evaluateN(-1,param,0,0,1,ncoeffs,sf_arr);
            return (qp+qm)/2.0 + constants[1]*(um-up)/(2.0);
            
            break;
    }
    
    cell = 0;
    return 0;
}


double Grid::R_u(int r,int c)
{
    Cell* cell = &grid[r][c];
    
    double volint = 0;
    double line_int1 = 0;
    double line_int2 = 0;
    double line_int3 = 0;
    double line_int4 = 0;
    double Ru = 0;
    double dt = cell->getDy();
    double dx = cell->getDx();
    
    //The volume residual term
    for(int i = 0; i < weights.length(); i++)
    {
        for(int j = 0; j < weights.length(); j++)
        {
            volint += weights[i]*weights[j]*
            fabs(cell->evaluateN(abscissa[i],abscissa[j],0,1,0,ncoeffs,sf_arr)+
                 constants[0]*cell->evaluateN(abscissa[i],abscissa[j],1,0,0,ncoeffs,sf_arr)+
                 cell->evaluateN(abscissa[i],abscissa[j],1,0,1,ncoeffs,sf_arr) );
        }
    }
    
    volint*=cell->getJacobian();
    
    //The line integral residual term
    
    //line integral along edge 1
    for(int i = 0; i < weights.length(); i++)
    {
        line_int1 += weights[i]*(fabs(u_star(abscissa[i],r,c,1) -
                                      cell->evaluateN(abscissa[i],-1,0,0,0,ncoeffs,sf_arr)) );
    }
    line_int1 *= .5*dx;
    
    //line integral along edge 2
    for(int i = 0; i < weights.length(); i++)
    {
        line_int2 += weights[i]*(constants[0]*fabs(u_star(abscissa[i],r,c,2) -
                                                   cell->evaluateN(1,abscissa[i],0,0,0,ncoeffs,sf_arr) )+
                                 fabs(q_star(abscissa[i],r,c,2) - 
                                      cell->evaluateN(1,abscissa[i],0,0,1,ncoeffs,sf_arr)));
    }
    line_int2 *= .5*dt;
    
    //line integral along edge 3
    for(int i = 0; i < weights.length(); i++)
    {
        line_int3 += weights[i]*(fabs(u_star(abscissa[i],r,c,3) - 
                                      cell->evaluateN(abscissa[i],1,0,0,0,ncoeffs,sf_arr)) );
    }
    line_int3 *= -.5*dx;
    
    //line integral along edge 4
    for(int i = 0; i < weights.length(); i++)
    {
        line_int4 += weights[i]*(constants[0]*fabs(u_star(abscissa[i],r,c,4) - 
                                                   cell->evaluateN(-1,abscissa[i],0,0,0,ncoeffs,sf_arr) )+
                                 fabs(q_star(abscissa[i],r,c,4) - 
                                      cell->evaluateN(-1,abscissa[i],0,0,1,ncoeffs,sf_arr)));
    }
    line_int4 *= -.5*dt;
    
    cell = 0;
    return volint+line_int1+line_int2+line_int3+line_int4;
}


double Grid::Denom_u(int r,int c)
{
    Cell* cell = &grid[r][c];
    
    double line_int1 = 0;
    double line_int2 = 0;
    double line_int3 = 0;
    double line_int4 = 0;
    double dt = cell->getDy();
    double dx = cell->getDx();
    
    //The line integral residual term
    
    //line integral along edge 1
    for(int i = 0; i < weights.length(); i++)
    {
        line_int1 += weights[i]*(fabs(u_star(abscissa[i],r,c,1) ) );
    }
    line_int1 *= dx*.5;
    
    //line integral along edge 2
    for(int i = 0; i < weights.length(); i++)
    {
        line_int2 += weights[i]*(constants[0]*fabs(u_star(abscissa[i],r,c,2) )+
                                 fabs(q_star(abscissa[i],r,c,2) ));
    }
    line_int2 *= dt*.5;
    
    //line integral along edge 3
    for(int i = 0; i < weights.length(); i++)
    {
        line_int3 += weights[i]*(fabs(u_star(abscissa[i],r,c,3) ) );
    }
    line_int3 *= -dx*.5;
    
    //line integral along edge 4
    for(int i = 0; i < weights.length(); i++)
    {
        line_int4 += weights[i]*(constants[0]*fabs(u_star(abscissa[i],r,c,4)  )+
                                 fabs(q_star(abscissa[i],r,c,4) ));
    }
    line_int4 *= -dt*.5;
    
    cell = 0;
    return line_int1+line_int2+line_int3+line_int4;
}


double Grid::R_q(int r,int c)
{
    Cell* cell = &grid[r][c];
    
    double volint = 0;
    double line_int1 = 0;
    double line_int2 = 0;
    double line_int3 = 0;
    double line_int4 = 0;
    double Rq = 0;
    double dt = cell->getDy();
    double dx = cell->getDx();
    
    //The volume residual term
    for(int i = 0; i < weights.length(); i++)
    {
        for(int j = 0; j < weights.length(); j++)
        {
            volint += weights[i]*weights[j]*
            fabs(cell->evaluateN(abscissa[i],abscissa[j],0,0,1,ncoeffs,sf_arr)+
                 constants[2]*cell->evaluateN(abscissa[i],abscissa[j],0,1,1,ncoeffs,sf_arr)+
                 constants[0]*constants[2]*cell->evaluateN(abscissa[i],abscissa[j],1,0,1,ncoeffs,sf_arr)+
                 constants[3]*cell->evaluateN(abscissa[i],abscissa[j],1,0,0,ncoeffs,sf_arr));
        }
    }
    
    volint*=cell->getJacobian();
    
    //The line integral residual term
    
    //line integral along edge 1
    for(int i = 0; i < weights.length(); i++)
    {
        line_int1 += weights[i]*constants[2]*(fabs(q_star(abscissa[i],r,c,1) - 
                                                   cell->evaluateN(abscissa[i],-1,0,0,1,ncoeffs,sf_arr)) );
    }
    line_int1 *= .5*dx;
    
    
    //line integral along edge 2
    for(int i = 0; i < weights.length(); i++)
    {
        line_int2 += weights[i]*(constants[3]*fabs(u_star(abscissa[i],r,c,2) - 
                                                   cell->evaluateN(1,abscissa[i],0,0,0,ncoeffs,sf_arr) )+
                                 constants[0]*constants[2]*fabs(q_star(abscissa[i],r,c,2) - 
                                                                cell->evaluateN(1,abscissa[i],0,0,1,ncoeffs,sf_arr)));
    }
    line_int2 *= .5*dt;
    
    
    //line integral along edge 3
    for(int i = 0; i < weights.length(); i++)
    {
        line_int3 += weights[i]*constants[2]*(fabs(q_star(abscissa[i],r,c,3) - 
                                                   cell->evaluateN(abscissa[i],1,0,0,1,ncoeffs,sf_arr)) );
    }
    line_int3 *= -.5*dx;
    
    //line integral along edge 4
    for(int i = 0; i < weights.length(); i++)
    {
        line_int4 += weights[i]*(constants[3]*fabs(u_star(abscissa[i],r,c,4) - 
                                                   cell->evaluateN(-1,abscissa[i],0,0,0,ncoeffs,sf_arr) )+
                                 constants[0]*constants[2]*fabs(q_star(abscissa[i],r,c,4) - 
                                                                cell->evaluateN(-1,abscissa[i],0,0,1,ncoeffs,sf_arr)));
    }
    line_int4 *= -.5*dt;
    
    
    cell = 0;
    return volint+line_int1+line_int2+line_int3+line_int4;
}


double Grid::Denom_q(int r,int c)
{
    
    Cell* cell = &grid[r][c];
    
    double dt = cell->getDy();
    double dx = cell->getDx();
    double line_int1 = 0;
    double line_int2 = 0;
    double line_int3 = 0;
    double line_int4 = 0;
    
    //The line integral residual term
    
    //line integral along edge 1
    for(int i = 0; i < weights.length(); i++)
    {
        line_int1 +=weights[i]*constants[2]*(fabs(q_star(abscissa[i],r,c,1)));
    }
    line_int1 *= .5*dx;
    
    
    //line integral along edge 2
    for(int i = 0; i < weights.length(); i++)
    {
        line_int2 += weights[i]*(constants[3]*fabs(u_star(abscissa[i],r,c,2) )+
                                 constants[0]*constants[2]*fabs(q_star(abscissa[i],r,c,2)) );
    }
    line_int2 *= .5*dt;
    
    
    //line integral along edge 3
    for(int i = 0; i < weights.length(); i++)
    {
        line_int3 += weights[i]*constants[2]*(fabs(q_star(abscissa[i],r,c,3) ) );
    }
    line_int3 *= -.5*dx;
    
    
    //line integral along edge 4
    for(int i = 0; i < weights.length(); i++)
    {
        line_int4 += weights[i]*(constants[3]*fabs(u_star(abscissa[i],r,c,4) )+
                                 constants[0]*constants[2]*fabs(q_star(abscissa[i],r,c,4) ));
    }
    line_int4 *= -.5*dt;
    
    
    cell = 0;
    return line_int1+line_int2+line_int3+line_int4;
}

double Grid::convergence_u()
{
    double Ru = 0;
    double Du = 0;
    
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            Ru += R_u(i,j);
            Du += Denom_u(i,j);
        }
    }
    
    return Ru/Du;
}

double Grid::convergence_q()
{
    double Rq = 0;
    double Dq = 0;
    
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            Rq += R_q(i,j);
            Dq += Denom_q(i,j);
        }
    }
    
    return Rq/Dq;
}