#pragma once
#include <Surface.h>
#include <map>
#include <Grid.h>
//It formulates the Poisson Disk Sampling for node generation according to following links:
//https://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf
//https://sighack.com/post/poisson-disk-sampling-bridsons-algorithm
class PoissonDiskSampling
{
private:
	Grid grid;

	void DiscretizeBoundaries(Surface*, double, double, int*, std::vector<pair<int, int>*>*, std::vector<pair<int, int>*>*, map<pair<int, int>, double*>*);
	bool IsValidPoint(Surface*, double*, double, double, int, int, std::vector<pair<int, int>*>*, std::vector<pair<int, int>*>*, map<pair<int, int>, double*>*);
	bool CheckCollisionWithBoundaries(std::vector<pair<int, int>*>*, double*, double, map<pair<int, int>, double*>*);
	bool CheckBoundaryDirection(map<pair<int, int>, double*>* , std::vector<std::pair<int,int>*>*, int, int);
	void ReorderBoundaryEdgesDirection(std::vector<std::pair<int, int>*>*, int, int);
	void LexicoGraphicNodeGeneration(map<pair<int, int>, double*>*, std::vector<std::pair<int, int>*>*, std::vector<std::pair<int, int>*>*);
public:
	PoissonDiskSampling(Surface*, double, int, int N = 2);
	Grid* PassGrid();
};

