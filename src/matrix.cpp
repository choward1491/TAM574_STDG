
#include "matrix.h"


/*
 * Helper function for clearing the memory 
 * of the matrix
 *
 *
 *
 */
void Matrix::clear()
{
	if( arr != NULL )
	{
		for(int i = 0; i < rows; i++)
			delete [] arr[i];
		delete [] arr;

		arr = NULL;
	}
	rows = 0; cols = 0;

}

Matrix::~Matrix()
{
	clear();
}
/*
 * Helper function for copying one matrix
 * to another
 *
 *
 *
 */
void Matrix::copy(const Matrix & m)
{
	arr = new double*[m.Rows()]();
	rows = m.Rows(); cols = m.Cols();
	for(int i = 0; i < rows; i++)
		arr[i] = new double[cols];

	for(int i = 0; i < rows; i ++)
	{
		for(int j = 0; j < cols; j++)
			arr[i][j] = m(i,j);
	}

}

Matrix::Matrix(const Matrix & m)
{
	copy(m);
}

/*
 * Create matrix based off of dimensions
 * 
 *
 *
 *
 */
Matrix::Matrix(int row, int col)
{
	rows = row; cols = col;

	arr = new double*[rows]();
	for( int i = 0; i < rows; i++)
		arr[i] = new double[cols]();
}

/*
 * Create matrix based off of dimensions
 * and value you want the whole matrix to have
 *
 *
 *
 */
Matrix::Matrix(int row, int col,double all_value)
{
	rows = row; cols = col;

	arr = new double*[rows]();
	for( int i = 0; i < rows; i++)
		arr[i] = new double[cols];

	for( int i = 0; i < rows; i++ )
	{
		for( int j = 0; j < cols; j++ )
			arr[i][j] = all_value;
	}
}

void Matrix::write2file(string filename)
{
	if( arr != 0 )
	{

		ofstream Fu;
		Fu.open(filename,ios::trunc);

		//Do file for u first
		for( int i = 0; i < rows; i++)
		{
				for( int j = 0; j < cols; j++ )
				{
					Fu<<arr[i][j]<<" ";
				}
				Fu<<"\n";
		}
		Fu.close();

	}
}


/*
 * Resize Matrix
 * 
 *
 *
 *
 */
void Matrix::resize(int r,int c)
{
	clear();
	rows = r; cols = c;

	arr = new double*[rows]();
	for( int i = 0; i < rows; i++)
		arr[i] = new double[cols];

	for( int i = 0; i < rows; i++ )
	{
		for( int j = 0; j < cols; j++ )
			arr[i][j] = 0;
	}

}


/*
 * Be able to add matrices together
 * 
 *
 *
 *
 */
Matrix Matrix::operator+(const Matrix & m) const
{
	if( rows == m.Rows() && cols == m.Cols() )
	{
		Matrix out(rows,cols);

		for(int i = 0; i < rows; i++ )
		{
			for(int j = 0; j < cols; j++)
				out(i,j) = arr[i][j] + m(i,j);
		}

		return out;
	}

	return Matrix();
}

/*
 * Be able to subtract matrices
 * 
 *
 *
 *
 */
Matrix Matrix::operator-(const Matrix & m) const
{
	if( rows == m.Rows() && cols == m.Cols() )
	{
		Matrix out(rows,cols);

		for(int i = 0; i < rows; i++ )
		{
			for(int j = 0; j < cols; j++)
				out(i,j) = arr[i][j] - m(i,j);
		}

		return out;
	}

	return Matrix();
}

/*
 * Set matrices equal to each other
 * 
 *
 *
 *
 */
Matrix & Matrix::operator=(const Matrix & m)
{
	if( &m != this)
	{
		clear();
		copy(m);
	}

	return *this;
}


/*
 * Multiply matrices together
 * 
 *
 *
 *
 */
Matrix Matrix::operator*(const Matrix & m) const
{
	if( cols == m.Rows() )
	{
		Matrix out(rows,m.Cols());

		for(int i = 0; i < rows; i++)
		{
			for( int j = 0; j < m.Cols(); j++)
			{
				for( int k = 0; k < cols; k++)
				{
					if( arr[i][k] != 0 && m(k,j) != 0 )
						out(i,j) = out(i,j) + arr[i][k]*m(k,j);
				}

			}
		}

		return out;
	}

	return Matrix();
}


/*
 * Print a matrix
 * 
 *
 *
 *
 */
void Matrix::print()
{
	if( arr != 0 )
	{
		for( int i = 0; i < rows; i++ )
		{   cout<<"|";
			for( int j = 0; j < cols; j++ )
				cout<<" "<<arr[i][j];
			cout<<" |\n";
		}
	}
}

/*
 * Access matrix contents based off of
 * its row and column coordinates
 *
 *
 *
 */
double Matrix::operator()(int r, int c) const
{
	return arr[r][c];
}

/*
 * Access value's location at row r and column c
 * so you can change it
 *
 *
 *
 */
double & Matrix::operator()(int r, int c)
{
	return arr[r][c];
}

/*
 * Use gaussian elimination to solve matrix
 * -Column matrix input and column vector output
 *
 *
 *
 */
Matrix Matrix::gaussSolve(Matrix b)
{
	return Matrix();
}


/*
 * Use gaussian elimination to solve matrix
 * -Vector input 
 *
 *
 *
 */
Matrix Matrix::gaussSolve(Vector b)
{
	return Matrix();
}


/*
 * Use gaussian elimination to solve matrix
 * -Vector input and output
 *
 *
 *
 */
Vector Matrix::gaussSolveV(Vector b)
{
	return Vector();
}


void Matrix::diag(int dim, double value)
{
	resize(dim,dim);
	for(int i = 0; i < dim; i++)
		arr[i][i] = value;
}


void Matrix::multrow(int row,double value)
{
	for(int i = 0; i < cols; i++)
		arr[row][i]*=value;
}

void Matrix::addrow(int row_toadd,double bycoeff,int row_addto)
{
	for(int i = 0; i < cols; i++)
		arr[row_addto][i] += arr[row_toadd][i]*bycoeff;
}

void Matrix::swaprow(int row1,int row2)
{
	double tmp = 0;
	for(int i = 0; i < cols; i++)
	{
		tmp = arr[row1][i];
		arr[row1][i] = arr[row2][i];
		arr[row2][i] = tmp;
	}
}
void Matrix::invert(Matrix & m)
{
	//print(); cout<<endl;

	double C = 0;
	if(rows == cols && rows >0 )
	{
		m.diag(rows,1);
		
		for(int i = 0; i < rows; i++)
		{
			//Do row swapping if necessary
			if( arr[i][i] == 0 )
			{
				for(int k = i+1; k<rows; k++)
				{
					if( arr[k][i] != 0 )
					{
						swaprow(k,i);
						m.swaprow(k,i);
						break;
					}
				}
			}

			//normalize the i,i component
			C = 1.0/arr[i][i];
			multrow(i,C);
			m.multrow(i,C);


			for(int j = 0; j < cols; j++)
			{
				if( j != i )
				{
					if( arr[j][i] != 0 )
					{
					C = -arr[j][i];
					addrow(i,C,j);
					m.addrow(i,C,j);
					}
				}
			}//End j loop

		}//End i loop
	}//End if statement
	//m.print();
}//End


Vector Matrix::operator*(const Vector & v) const
{
	Vector output;
	if( cols == v.length() )
	{
		output.change(cols);

		for(int i = 0; i < rows; i++)
		{
			for(int j = 0; j < cols; j++)
				output[i]+=arr[i][j]*v[j];
		}
	}

	return output;
}












Upper::Upper(int d)
{
	dim = d;

	arr = new double*[dim];

	for(int i = 0; i < dim; i++)
	{
		arr[i] = new double[dim-i]();
	}
}

Upper::Upper(int d,double value)
{
	dim = d;

	arr = new double*[dim];

	for(int i = 0; i < dim; i++)
	{
		arr[i] = new double[dim-i]();
	}

	for(int i = 0; i < dim; i++)
	{
		for(int j = 0; j < dim-i; j++ )
			arr[i][j] = value;
	}
}

void Upper::clear()
{
	if( arr != 0 )
	{
	for(int i = 0; i < dim; i++)
		delete [] arr[i];
	delete [] arr;
	arr = 0;
	}
	dim = 0;

}
Upper::~Upper()
{
	clear();
}

void Upper::resize(int d)
{

	clear();
	dim = d;

	arr = new double*[dim];

	for(int i = 0; i < dim; i++)
	{
		arr[i] = new double[dim-i]();
	}

	for(int i = 0; i < dim; i++)
	{
		for(int j = 0; j < dim-i; j++ )
			arr[i][j] = 0;
	}
}

void Upper::print()
{
	if( arr != 0 )
	{
		for(int i = 0; i < dim; i++)
		{
			cout<<"| ";
			for(int j = 0; j < dim; j++ )
			{
				if( j >= i )
					cout<<arr[i][j-i]<<" ";
				else
					cout<<0<<" ";
			}
			cout<<"|\n";
		}

	}
}

void Upper::toDiagonal(double value)
{

	if( arr != NULL )
	{
		for(int i = 0; i < dim; i++)
		{
			for(int j = 0; j < dim-i; j++ ){

				if( j == 0 )
					arr[i][j] = value;
				else
					arr[i][j] = 0;
			}//j loop
		}//i loop
	}

}
double Upper::determinant()
{
	double output = 1;

	for(int i = 0; i < dim; i++)
		output *= arr[i][0];
	return output;
}

Vector Upper::solve(Vector const & b )
{
	Vector output;

	if( b.length() == dim && determinant() != 0)
	{
		output.change(dim);
		double sum = 0;

		for(int i = dim-1; i >= 0; i-- )
		{
			for(int j = dim-1-i,count = dim-1; j > 0; j--,count-- )
				sum += arr[i][j]*output[count];
			output[i] = (b[i]-sum)/arr[i][0];
			sum = 0;
		}
	}
	return output;
}

double Upper::operator()(int r, int c) const
{
	if( r > dim || c > dim || c > r || r < 0 || c < 0)
		return 0;
	else
		return arr[r][c];
}

double & Upper::operator()(int r, int c)
{
	if( r > dim || c > dim || c < r || r < 0 || c < 0){
		cout<<"Trying to access upper matrix element outside of bounds\n";
		return arr[0][0];
	}else
		return arr[r][c-r];
}







Lower::Lower(int d)
{
	dim = d;

	arr = new double*[dim];

	for(int i = 0; i < dim; i++)
	{
		arr[i] = new double[i+1]();
	}
}

Lower::Lower(int d,double value)
{
	dim = d;

	arr = new double*[dim];

	for(int i = 0; i < dim; i++)
	{
		arr[i] = new double[i+1]();
	}

	for(int i = 0; i < dim; i++)
	{
		for(int j = 0; j <=i; j++)
			arr[i][j] = value;
	}
}

void Lower::clear()
{
	if( arr != 0 ){
	for(int i = 0; i < dim; i++)
		delete [] arr[i];
	delete [] arr;
	}
	dim = 0;
}

Lower::~Lower()
{
	clear();
}

void Lower::print()
{
	if( arr != 0 )
	{
		for(int i = 0; i < dim; i++)
		{
			cout<<"| ";
			for(int j = 0; j < dim; j++ )
			{
				if( j <= i )
					cout<<arr[i][j]<<" ";
				else
					cout<<0<<" ";
			}
			cout<<"|\n";
		}

	}
}

void Lower::toDiagonal(double val)
{
	if( arr != 0 )
	{
			for(int i = 0; i < dim; i++)
			{
				for(int j = 0; j <=i; j++){
					if( i == j)
						arr[i][j] = val;
					else
						arr[i][j] = 0;
				}
			}
	}
}

double Lower::determinant()
{
	double output = 1;
	if( arr != 0 )
	{
		for(int i = 0; i < dim; i++)
			output*=arr[i][i];
		return output;
	}else
		return 0;
}

Vector Lower::solve(Vector const & b)
{
	Vector output;

	if( b.length() == dim && determinant() != 0)
	{
		output.change(dim);
		double sum = 0;

		for(int i = 0; i < dim; i++ )
		{
			for(int j = 0; j < i; j++ )
				sum += arr[i][j]*output[j];
			output[i] = (b[i]-sum)/arr[i][i];
			sum = 0;
		}
	}
	return output;
}

double Lower::operator()(int r, int c) const
{
	if( r > dim || c > dim || c > r || r < 0 || c < 0)
		return 0;
	else
		return arr[r][c];
}

double & Lower::operator()(int r, int c)
{
	if( r > dim || c > dim || c > r || r < 0 || c < 0){
		cout<<"Trying to access lower matrix element outside of bounds\n";
		return arr[0][0];
	}else
		return arr[r][c];
}

void Lower::resize(int d)
{
	clear();

	dim = d;

	arr = new double*[dim];

	for(int i = 0; i < dim; i++)
	{
		arr[i] = new double[i+1]();
	}

	for(int i = 0; i < dim; i++)
	{
		for(int j = 0; j <=i; j++)
			arr[i][j] = 0;
	}
}


void LU_Decomposition(Lower & L, Upper & U, const Matrix & A)
{
	if( A.Cols() != A.Rows() ){
		cout<<"The matrix isn't square, so no LU decomposition available\n";
		return;
	}

	int dim = A.Cols();

	L.resize(dim); U.resize(dim);
	L.toDiagonal(1);
	double sum = 0;

	for(int i = 0; i < dim; i++)
	{
		for(int j = 0; j < dim; j++ )
		{
			if( j < i )
			{
				if( U(j,j) == 0 )
				{
					L(i,j) = 1;
				}else{
					for( int k = 0; k <= j-1; k++)
						sum += L(i,k)*U(k,j);
					L(i,j) = (A(i,j)-sum)/U(j,j); sum = 0;
				}

			}else{
				for( int k = 0; k <= i-1; k++ )
					sum += L(i,k)*U(k,j);
				U(i,j) = A(i,j) - sum; sum = 0;
			}//End of if statement
		}//End of j for loop
	}//End of i for loop

	//End
}