#pragma once
#include "IDelaunayTriangulation.h"
#include <map>
#include <vector>
class DelaunayTriangulation2D : public IDelaunayTriangulation
{
private:

	//std::vector<std::pair<int, int>> segments;
	//std::vector<int*> elements;
	//std::map<int, std::pair<double, double>>* vertices{ nullptr };

	int* FinishEdge(std::pair<int, int>*);
	bool CheckSegmentsCollision(double*, double*, double*, double*);
	void AddTriangle(int*, int*);
	void PassSegments(std::vector<std::pair<int, int>>*, std::vector<std::pair<int, int>>*);
	void PassVertices(std::map<int, std::pair<double, double>>*);
protected:
	virtual double Orient(double*) override final;
	virtual bool InsideArea(double*) override final;
public:
	//Grid grid;

	DelaunayTriangulation2D(Grid*/*, std::vector<std::pair<int, int>>*, std::vector<std::pair<int, int>>**/);
};

