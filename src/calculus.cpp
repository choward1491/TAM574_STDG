
#include "calculus.h"

int simplesum(int start, int end,int increment)
{
	int tmp = 0;
	for( int i = start; i <=end; i+=increment)
		tmp = tmp+i;
	return tmp;
}

Vector getQuadrature(int p)
{
   
   int cnt = 1;
   int n = ceil( (2*p+1) / 2.0 );
   Vector out(2*n);
   string line;

   //Get the Abscissa
   ifstream myfile ("GaussAbscissa.txt");
   if (myfile.is_open())
   {
    while ( myfile.good() )
    {
	   while( cnt < n && getline( myfile, line )  )  // read each line:
	   {
		   cnt++;
	   }

	   getline( myfile,line);
	   istringstream is( line );
	   double val; int cnt1 = 0;
	   while( is >> val ) {   // read each number in line
		   //cout<<val<<endl;
		   out[cnt1] = val; cnt1++;
	   }

	   break;
    }
    myfile.close();
   }
    cnt = 1;
   //Get the weights
   ifstream myfile2 ("GaussWeights.txt");
   if (myfile2.is_open())
   {
    while ( myfile2.good() )
    {
	   while( cnt != n && getline( myfile2, line ) )  // read each line:
	   {
		   cnt++;
	   }

	   getline( myfile2,line);
	   istringstream is( line );
	   double val; int cnt1 = 0;
	   while( is >> val ) {   // read each number in line
		   out[cnt1+n] = val; cnt1++;
	   }

	   break;
    }
    myfile2.close();
   }

   return out;
}
