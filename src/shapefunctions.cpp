
#include "shapefunctions.h"

Vector getNormedShapeCoeffs(int p,double dx, double dt)
{
	int tmp = simplesum(1,p+1,1);

	Vector out(tmp);
	int cnt = 1;
	int curr = 0;
	int alpha = 0;
	int beta = 0;
	for(int i = 0; i < tmp; i++)
	{
		beta = i-curr;
		alpha = (cnt-1)-beta;

		//cout<<"alpha = "<<alpha<<", beta = "<<beta<<"\n";

		if( i == curr+cnt-1){
			curr+=cnt;
			cnt++;
		}

		out[i] = sqrt( (2*alpha+1)*(2*beta+1)/(dx*dt) );
	}

	return out;
}

Points getShapeFuncPolyOrders(int p){

	int tmp = simplesum(1,p+1,1);

	Points out(tmp);
	int cnt = 1;
	int curr = 0;
	int alpha = 0;
	int beta = 0;
	for(int i = 0; i < tmp; i++)
	{
		beta = i-curr;
		alpha = (cnt-1)-beta;

		if( i == curr+cnt-1){
			curr+=cnt;
			cnt++;
		}

		out[i].x = alpha; out[i].y = beta;
	}

	return out;
}