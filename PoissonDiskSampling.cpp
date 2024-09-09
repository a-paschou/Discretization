#include "PoissonDiskSampling.h"
#include<cstdlib>
#include <time.h>
#include <iostream>
//TODO: Find a better way to determine whether a new node is located between the boundaries
PoissonDiskSampling::PoissonDiskSampling(Surface* surface, double radius, int nAttempts, int N)
{
	double cellSize{ radius };
	radius *= sqrt(N);
	double width{ surface->GetMaxX() - surface->GetMinX()};
	double height{ surface->GetMaxY() - surface->GetMinY()};
	int nCellsWidth{ (int)ceil(width / cellSize) + 1 };
	int nCellsHeight{ (int)ceil(height / cellSize) + 1 };
	std::vector<pair<int, int>*> auxOuterBoundaryEdges;
	std::vector<pair<int, int>*> auxInnerBoundaryEdges;
	map<pair<int, int>, double*> auxNodesMap;
	std::vector<double*> activeNodes;

	int nodeID{ 0 };
	DiscretizeBoundaries(surface, radius, cellSize, &nodeID, &auxOuterBoundaryEdges, &auxInnerBoundaryEdges, &auxNodesMap);

	srand(time(0));
	auto it{ auxNodesMap.begin() };
	std::advance(it, rand() % auxNodesMap.size());
	activeNodes.push_back(new double[2] {it->second[0], it->second[1]});
	double* pointCoors{ nullptr };
	pointCoors = new double[3];
	const double pi{ 3.14159265358979323846 };
	double newRadius{ 0 };
	double theta{ 0 };
	bool found{};

	while (activeNodes.size() > 0)
	{
		found = false;
		int randomIndex{ (int)(rand() % (activeNodes.size())) };
		pointCoors[0] = activeNodes.at(randomIndex)[0];
		pointCoors[1] = activeNodes.at(randomIndex)[1];
		for (size_t i = 0; i < nAttempts; i++)
		{
			theta = ((double)rand() / RAND_MAX) * (2 * pi);
			newRadius = radius + ((double)rand() / RAND_MAX) * (2 * radius - radius);
			double newCoor[2]{ pointCoors[0] + newRadius * cos(theta), pointCoors[1] + newRadius * sin(theta) };
			if (!IsValidPoint(surface, newCoor, cellSize, radius, nCellsWidth, nCellsHeight, &auxOuterBoundaryEdges, &auxInnerBoundaryEdges, &auxNodesMap))
			{
				continue;
			}
			auxNodesMap[make_pair((int)(floor(newCoor[0] / cellSize)), (int)(floor(newCoor[1] / cellSize)))] = new double[2] {newCoor[0], newCoor[1]};
			vertices[nodeID] = make_pair(newCoor[0], newCoor[1]);
			grid.nodesMap[nodeID] = new Node(nodeID, newCoor[0], newCoor[1]);
			nodeID++;
			activeNodes.push_back(new double[2] { newCoor[0], newCoor[1] });
			found = true;
			break;
		}
		if (!found)
		{
			delete activeNodes.at(randomIndex);
			activeNodes.erase(activeNodes.begin() + randomIndex);
		}

	}
	delete[] pointCoors;


	//auto it2{ vertices.begin() };
	//for (it2 = vertices.begin(); it2 != vertices.end(); it2++)
	//{
	//	std::cout << it2->second.first << " " << it2->second.second << std::endl;
	//}


	auxNodesMap.clear();
	auxOuterBoundaryEdges.clear();
	auxInnerBoundaryEdges.clear();
}

Grid* PoissonDiskSampling::PassGrid()
{
	return &grid;
}

void PoissonDiskSampling::DiscretizeBoundaries(Surface* surface, double radius, double cellSize, int* nodeID, std::vector<pair<int, int>*>* auxOuterBoundaryEdges, std::vector<pair<int, int>*>* auxInnerBoundaryEdges, map<pair<int, int>, double*>* auxNodesMap)
{
	int currentBoundary{ 0 };
	int numberOfBoundaries{ (int)surface->GetOuterBoundary()->size() };
	int initialNodeID{};
	pair<int, int> currentPosition{};
	pair<int, int> previousPosition{};
	int currentNodeID{ *nodeID };
	int previousNodeID{};
	pair<int, int> initialPosition{};
	double* coor{ nullptr };
	coor = new double[2];
	for (IShape* boundary : *surface->GetOuterBoundary())
	{
		int nEdges{ (int)(floor(boundary->CalculateLength() / radius)) };
		double newRadius{ boundary->CalculateLength() / nEdges };
		for (size_t i = 0; i <= nEdges; i++)
		{
			if (i == 0)
			{
				if (currentBoundary == 0)
				{
					coor = boundary->GetStartPoint();
					currentPosition = make_pair(floor(coor[0] / cellSize), floor(coor[1] / cellSize));
					//outerBoundaryNodes.push_back(&currentPosition);
					(*auxNodesMap)[currentPosition] = new double[2] { coor[0], coor[1] };
					vertices[*nodeID] = make_pair(coor[0], coor[1]);
					grid.nodesMap[*nodeID] = new Node(*nodeID, coor[0], coor[1]);
					previousNodeID = *nodeID;
					initialNodeID = *nodeID;
					*nodeID += 1;
					previousPosition = currentPosition;
					initialPosition = previousPosition;
				}
			}
			else if (i < nEdges)
			{
				coor = boundary->CalculatePointOnDistance(newRadius * i);
				currentPosition = make_pair(floor(coor[0] / cellSize), floor(coor[1] / cellSize));

				(*auxNodesMap)[currentPosition] = new double[2] { coor[0], coor[1] };
				currentNodeID = *nodeID;
				vertices[currentNodeID] = make_pair(coor[0], coor[1]);
				grid.nodesMap[*nodeID] = new Node(*nodeID, coor[0], coor[1]);
				*nodeID += 1;

				pair<int, int>* newPair = new pair<int, int>[2] {previousPosition, currentPosition};
				//outerBoundaryEdges.push_back(make_pair(previousNodeID, currentNodeID));
				grid.segments.push_back(make_pair(previousNodeID, currentNodeID));
				previousNodeID = currentNodeID;
				auxOuterBoundaryEdges->push_back(newPair);
				previousPosition = currentPosition;
			}
			else
			{
				if (currentBoundary < numberOfBoundaries - 1)
				{
					coor = boundary->GetEndPoint();
					currentPosition = make_pair(floor(coor[0] / cellSize), floor(coor[1] / cellSize));

					(*auxNodesMap)[currentPosition] = new double[2] { coor[0], coor[1] };
					currentNodeID = *nodeID;
					vertices[currentNodeID] = make_pair(coor[0], coor[1]);
					grid.nodesMap[*nodeID] = new Node(*nodeID, coor[0], coor[1]);
					*nodeID += 1;

					pair<int, int>* newPair = new pair<int, int>[2] {previousPosition, currentPosition};
					//outerBoundaryEdges.push_back(make_pair(previousNodeID, currentNodeID));
					grid.segments.push_back(make_pair(previousNodeID, currentNodeID));
					auxOuterBoundaryEdges->push_back(newPair);
					previousNodeID = currentNodeID;
					previousPosition = currentPosition;
				}
				else
				{
					//outerBoundaryEdges.push_back(make_pair(previousNodeID, initialNodeID));
					grid.segments.push_back(make_pair(previousNodeID, initialNodeID));
					pair<int, int>* newPair = new pair<int, int>[2] {previousPosition, initialPosition};
					auxOuterBoundaryEdges->push_back(newPair);
				}
			}
		}
		currentBoundary++;
	}
	if (CheckBoundaryDirection(0, grid.segments.size() - 1))
	{
		ReorderBoundaryEdgesDirection(0, outerBoundaryEdges.size() - 1);
	}

	int firstBoundariesEdge{};
	for (std::vector<IShape*>* boundaries : *surface->GetInnerBoundary())
	{
		currentBoundary = 0;
		numberOfBoundaries = boundaries->size();
		pair<int, int> currentPosition{};
		pair<int, int> previousPosition{};
		int currentNodeID{ *nodeID };
		int previousNodeID{};
		firstBoundariesEdge = grid.segments.size();
		for (IShape* boundary : *boundaries)
		{
			int nEdges{ (int)(floor(boundary->CalculateLength() / radius)) };
			double newRadius{ boundary->CalculateLength() / nEdges };
			for (size_t i = 0; i <= nEdges; i++)
			{
				if (i == 0)
				{
					if (currentBoundary == 0)
					{
						coor = boundary->GetStartPoint();
						currentPosition = make_pair(floor(coor[0] / cellSize), floor(coor[1] / cellSize));

						(*auxNodesMap)[currentPosition] = new double[2] { coor[0], coor[1] };;
						vertices[*nodeID] = make_pair(coor[0], coor[1]);
						grid.nodesMap[*nodeID] = new Node(*nodeID, coor[0], coor[1]);
						previousNodeID = *nodeID;
						initialNodeID = *nodeID;
						*nodeID += 1;
						previousPosition = currentPosition;
						initialPosition = previousPosition;
					}
				}
				else if (i < nEdges)
				{
					coor = boundary->CalculatePointOnDistance(newRadius * i);
					currentPosition = make_pair(floor(coor[0] / cellSize), floor(coor[1] / cellSize));

					(*auxNodesMap)[currentPosition] = new double[2] { coor[0], coor[1] };
					currentNodeID = *nodeID;
					vertices[currentNodeID] = make_pair(coor[0], coor[1]);
					grid.nodesMap[*nodeID] = new Node(*nodeID, coor[0], coor[1]);
					*nodeID += 1;

					pair<int, int>* newPair = new pair<int, int>[2] {previousPosition, currentPosition};
					//innerBoundaryEdges.push_back(make_pair(previousNodeID, currentNodeID));
					grid.segments.push_back(make_pair(previousNodeID, currentNodeID));
					previousNodeID = currentNodeID;
					auxInnerBoundaryEdges->push_back(newPair);
					previousPosition = currentPosition;
				}
				else
				{
					if (currentBoundary < numberOfBoundaries - 1)
					{
						coor = boundary->GetEndPoint();
						currentPosition = make_pair(floor(coor[0] / cellSize), floor(coor[1] / cellSize));
						(*auxNodesMap)[currentPosition] = new double[2] { coor[0], coor[1] };;
						currentNodeID = *nodeID;
						vertices[currentNodeID] = make_pair(coor[0], coor[1]);
						grid.nodesMap[*nodeID] = new Node(*nodeID, coor[0], coor[1]);
						*nodeID += 1;

						pair<int, int>* newPair = new pair<int, int>[2] {previousPosition, currentPosition};
						//innerBoundaryEdges.push_back(make_pair(previousNodeID, currentNodeID));
						grid.segments.push_back(make_pair(previousNodeID, currentNodeID));
						auxInnerBoundaryEdges->push_back(newPair);
						previousNodeID = currentNodeID;
						previousPosition = currentPosition;
					}
					else
					{
						//innerBoundaryEdges.push_back(make_pair(previousNodeID, initialNodeID));
						grid.segments.push_back(make_pair(previousNodeID, initialNodeID));
						pair<int, int>* newPair = new pair<int, int>[2] {previousPosition, initialPosition};
						auxInnerBoundaryEdges->push_back(newPair);
					}
				}
			}
			currentBoundary++;
		}
		if (!CheckBoundaryDirection(firstBoundariesEdge, grid.segments.size() - 1))
		{
			ReorderBoundaryEdgesDirection(firstBoundariesEdge, grid.segments.size() - 1);
		}
	}
	delete[] coor;
}

bool PoissonDiskSampling::IsValidPoint(Surface* surface, double* coor, double cellSize, double radius, int nCellWidth, int nCellHeight, std::vector<pair<int, int>*>* auxOuterBoundaryEdges, std::vector<pair<int, int>*>* auxInnerBoundaryEdges, map<pair<int, int>, double*>* auxNodesMap)
{
	if (!CheckCollisionWithBoundaries(auxOuterBoundaryEdges, coor, (surface->GetMinX() - 1), auxNodesMap) || CheckCollisionWithBoundaries(auxInnerBoundaryEdges, coor, (surface->GetMinX() - 1), auxNodesMap))
	{
		return false;
	}
	int xIndex{ (int)floor(coor[0] / cellSize) };
	int yIndex{ (int)floor(coor[1] / cellSize) };
	int i0{ max(xIndex - 1, 0) };
	int i1{ min(xIndex + 1, nCellWidth - 1) };
	int j0{ max(yIndex - 1, 0) };
	int j1{ min(yIndex + 1, nCellHeight - 1) };

	for (int i = i0; i <= i1; i++)
	{
		for (int j = j0; j <= j1; j++)
		{
			if (auxNodesMap->find(make_pair(i, j)) != auxNodesMap->end())
			{
				double distance{ sqrt(pow(auxNodesMap->at(make_pair(i, j))[0] - coor[0], 2) + pow(auxNodesMap->at(make_pair(i, j))[1] - coor[1], 2)) };
				if (distance < radius)
				{
					return false;
				}
			}
		}
	}
	return true;
}

bool PoissonDiskSampling::CheckCollisionWithBoundaries(std::vector<pair<int, int>*>* edges, double* coor, double auxiliaryCoor, map<pair<int, int>, double*>* nodesMap)
{
	int numberOfIntersects{ 0 };
	double* const auxiliaryPoint{ new double[2] };
	auxiliaryPoint[0] = auxiliaryCoor;
	auxiliaryPoint[1] = coor[1];
	double* const pointStart{ new double[2] };
	//pointStart = new double[2];
	double* const pointEnd{ new double[2] };
	//pointEnd = new double[2];
	for (size_t i = 0; i < edges->size(); i++)
	{
		pointStart[0] = nodesMap->at(edges->at(i)[0])[0];
		pointStart[1] = nodesMap->at(edges->at(i)[0])[1];
		pointEnd[0] = nodesMap->at(edges->at(i)[1])[0];
		pointEnd[1] = nodesMap->at(edges->at(i)[1])[1];

		double beta2{ coor[1] };
		double intersectPoint[2]{ 0,0 };

		double minX{};
		double maxX{};
		double minY{};
		double maxY{};
		if ((pointStart[0] != pointEnd[0]))
		{
			if ((abs(pointStart[1] - pointEnd[1]) >= pow(10, -10)))
			{
				double alpha1{ (pointStart[1] - pointEnd[1]) / (pointStart[0] - pointEnd[0]) };
				double beta1{ pointStart[1] - alpha1 * pointStart[0] };
				intersectPoint[0] = (beta2 - beta1) / (alpha1 - 0);
				intersectPoint[1] = alpha1 * intersectPoint[0] + beta1;
				minX = min(coor[0], auxiliaryPoint[0]);
				maxX = max(coor[0], auxiliaryPoint[0]);
				minY = min(pointStart[1], pointEnd[1]);
				maxY = max(pointStart[1], pointEnd[1]);
			}
			else
			{
				continue;
			}
		}
		else
		{
			intersectPoint[0] = pointStart[0];
			intersectPoint[1] = beta2;
			minX = min(coor[0], auxiliaryPoint[0]);
			maxX = max(coor[0], auxiliaryPoint[0]);
			minY = min(pointStart[1], pointEnd[1]);
			maxY = max(pointStart[1], pointEnd[1]);
		}

		if (intersectPoint[0] >= minX && intersectPoint[0] <= maxX && intersectPoint[1] >= minY && intersectPoint[1] <= maxY)
		{
			numberOfIntersects++;
		}
	}
	delete[] pointStart;
	delete[] pointEnd;
	delete[] auxiliaryPoint;
	int result{ numberOfIntersects % 2 };
	if (result == 0)
	{
		return false;
	}
	else
	{
		return true;
	}
}

// Determines whether the points are in clockwise (true) or counter-clockwise (false) order.
bool PoissonDiskSampling::CheckBoundaryDirection(int firstEdge, int lastEdge)
{
	double sum{ 0 };
	for (size_t i = firstEdge; i <= lastEdge; i++)
	{
		sum += (vertices.at(grid.segments.at(i).second).first - vertices.at(grid.segments.at(i).first).first) * (vertices.at(grid.segments.at(i).second).second + vertices.at(grid.segments.at(i).first).second);
	}
	if (sum > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

// Interchanges the boundary's edges first and second vetice
void PoissonDiskSampling::ReorderBoundaryEdgesDirection(int firstEdge, int lastEdge)
{
	int temp{};
	for (size_t i = firstEdge; i <= lastEdge; i++)
	{
		temp = grid.segments.at(i).first;
		//boundary->at(i).first = boundary->at(i).second;
		//boundary->at(i).second = temp;
		grid.segments.at(i).first = grid.segments.at(i).second;
		grid.segments.at(i).second = temp;
	}
}
