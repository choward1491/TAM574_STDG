#include "geometry.h"

Points::Points(int num)
{
	count = num;
	pts = new Point[num]();
}

void Points::clear()
{
	if( pts == 0 )
		delete [] pts;
	pts = 0;
	count = 0;
}

Points::~Points()
{
	clear();
}

void Points::changePoints(int num)
{
	clear();
	count = num;
	pts = new Point[num]();
}