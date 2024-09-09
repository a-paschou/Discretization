#pragma once
#include <string>
#include "Grid.h"
class IDelaunayTriangulation
{
protected:
	//Checks whether a number of points are arranged clockwise, counter-clockwise or they are colinear
	virtual double Orient(double*) = 0;
	//Checks whether a point lies inside a circumcircle (2D) or a circumscribed sphere (3D)
	virtual bool InsideArea(double*) = 0;
public:
	//virtual std::string ExportGrid() = 0;
	Grid* grid;
};

