


#ifndef GEOMETRY_H
#define GEOMETRY_H


#include <iostream>
using namespace std;

class Point
{
	
public:
	Point():x(0),y(0){}
	Point(double X, double Y):x(X),y(Y){}
	void change(double X, double Y){x = X; y = Y;}
	void print(){ cout<<"( "<<x<<", "<<y<<") "; }
	double x;
	double y;
};


class Points
{
public:
	Points():count(0),pts(0){}
	Points(int num);
	~Points();
	void changePoints(int num);
	Point operator[](int index) const { return pts[index];}
	Point & operator[](int index){ return pts[index];}
	void print()
	{   cout<<"[ ";
		for(int i = 0; i < count; i++){
			pts[i].print(); cout<<", ";
		}   cout<<" ]\n";
	}

private:
	Point* pts;
	int count;
	void clear();

};
#endif