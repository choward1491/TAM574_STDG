

#ifndef SHAPEFUNCTIONS_H
#define SHAPEFUNCTIONS_H

#include "vector.h"
#include "calculus.h"
#include "geometry.h"

Vector getNormedShapeCoeffs(int p,double dx, double dt);
Points getShapeFuncPolyOrders(int p);

#endif