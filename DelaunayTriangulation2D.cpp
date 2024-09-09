#include "DelaunayTriangulation2D.h"
#include <corecrt_math.h>
#include <iostream>
#include "CST.h"

DelaunayTriangulation2D::DelaunayTriangulation2D(Grid* grid/*, std::vector<std::pair<int, int>>* outerBoundaryEdges, std::vector<std::pair<int, int>>* innerBoundaryEdges*/)
{
    //grid.vertices = vertices;
    this->grid = grid;
    //PassVertices(vertices);
    //PassSegments(outerBoundaryEdges, innerBoundaryEdges);
    int elementID{ 0 };
    std::pair<int, int>* edge{ nullptr };
    std::pair<int, int>* edge2{ nullptr };
    std::pair<int, int>* edge3{ nullptr };
    edge = new std::pair<int, int>();
    edge2 = new std::pair<int, int>();
    edge3 = new std::pair<int, int>();
    while (this->grid->segments.size() > 0)
    {
        edge[0] = std::make_pair(this->grid->segments[0].first, this->grid->segments[0].second);
        this->grid->segments.erase(grid->segments.begin());
        int* triangleVertices{ nullptr };
        triangleVertices = new int[3];
        triangleVertices = FinishEdge(edge);
        if (triangleVertices != nullptr)
        {
            AddTriangle(&elementID, triangleVertices);
            edge2->first = this->grid->elementsMap.at(elementID - 1)->GetNodeID(2);
            edge2->second = this->grid->elementsMap.at(elementID - 1)->GetNodeID(0);
            edge3->first = this->grid->elementsMap.at(elementID - 1)->GetNodeID(1);
            edge3->second = this->grid->elementsMap.at(elementID - 1)->GetNodeID(2);

            //edge2->first = grid.elements.at(elementID - 1)[2];
            //edge2->second = grid.elements.at(elementID - 1)[0];
            //edge3->first = grid.elements.at(elementID - 1)[1];
            //edge3->second = grid.elements.at(elementID - 1)[2];
            bool flagEdge2{ false };
            bool flagEdge3{ false };
            int positionEdge2{ 0 };
            int positionEdge3{ 0 };
            for (size_t i = 0; i < this->grid->segments.size(); i++)
            {
                if (!flagEdge2)
                {
                    if ((this->grid->segments.at(i).first == edge2->first && this->grid->segments.at(i).second == edge2->second) || (this->grid->segments.at(i).first == edge2->second && this->grid->segments.at(i).second == edge2->first))
                    {
                        flagEdge2 = true;
                        positionEdge2 = i;
                    }
                }
                if (!flagEdge3)
                {
                    if ((this->grid->segments.at(i).first == edge3->first && this->grid->segments.at(i).second == edge3->second) || (this->grid->segments.at(i).first == edge3->second && this->grid->segments.at(i).second == edge3->first))
                    {
                        flagEdge3 = true;
                        positionEdge3 = i;
                    }
                }
                if (flagEdge2 && flagEdge3)
                {
                    break;
                }
            }
            if (flagEdge2 && flagEdge3)
            {
                this->grid->segments.erase(this->grid->segments.begin() + std::max(positionEdge2, positionEdge3));
                this->grid->segments.erase(this->grid->segments.begin() + std::min(positionEdge2, positionEdge3));
            }
            else if (flagEdge2)
            {
                this->grid->segments.erase(this->grid->segments.begin() + positionEdge2);
            }
            else if (flagEdge3)
            {
                this->grid->segments.erase(this->grid->segments.begin() + positionEdge3);
            }
            if (!flagEdge2)
            {
                int temp1{ edge2->first };
                int temp2{ edge2->second };
                edge2->first = temp2;
                edge2->second = temp1;
                this->grid->segments.push_back(*edge2);
            }
            if (!flagEdge3)
            {
                int temp1{ edge3->first };
                int temp2{ edge3->second };
                edge3->first = temp2;
                edge3->second = temp1;
                this->grid->segments.push_back(*edge3);
            }
        }
    }
}

int* DelaunayTriangulation2D::FinishEdge(std::pair<int, int>* edge)
{
    int points[3]{};
    points[0] = edge->first;
    points[1] = edge->second;

    double* coors{ nullptr };
    coors = new double[8];

    coors[0] = grid->nodesMap.at(points[0])->GetX();
    coors[1] = grid->nodesMap.at(points[0])->GetY();
    coors[2] = grid->nodesMap.at(points[1])->GetX();
    coors[3] = grid->nodesMap.at(points[1])->GetY();


    //coors[0] = grid.vertices->at(points[0]).first;
    //coors[1] = grid.vertices->at(points[0]).second;
    //coors[2] = grid.vertices->at(points[1]).first;
    //coors[3] = grid.vertices->at(points[1]).second;
    double* tryCoors{ nullptr };
    tryCoors = new double[6];
    tryCoors[0] = grid->nodesMap.at(points[0])->GetX();
    tryCoors[1] = grid->nodesMap.at(points[0])->GetY();
    tryCoors[2] = grid->nodesMap.at(points[1])->GetX();
    tryCoors[3] = grid->nodesMap.at(points[1])->GetY();

    //tryCoors[0] = grid.vertices->at(points[0]).first;
    //tryCoors[1] = grid.vertices->at(points[0]).second;
    //tryCoors[2] = grid.vertices->at(points[1]).first;
    //tryCoors[3] = grid.vertices->at(points[1]).second;
    double midPoint[2]{ (tryCoors[0] + tryCoors[2]) / 2, (tryCoors[1] + tryCoors[3]) / 2 };

    bool firstTriangleCreated{ false };

    auto it{ grid->nodesMap.begin() };
    for (it = grid->nodesMap.begin(); it != grid->nodesMap.end(); it++)
    {
        if (it->first == points[0] || it->first == points[1])
        {
            continue;
        }
        tryCoors[4] = it->second->GetX();
        tryCoors[5] = it->second->GetY();
        if (Orient(tryCoors) > 0)
        {
            coors[6] = tryCoors[4];
            coors[7] = tryCoors[5];
            if (!firstTriangleCreated || InsideArea(coors))
            {
                bool flag{ true };
                for (auto &edgeVertices: grid->segments)
                {
                    if (it->first == edgeVertices.first || it->first == edgeVertices.second)
                    {
                        continue;
                    }
                    double segmentPointStart[2]{ 0,0 };
                    double segmentPointEnd[2]{ 0,0 };
                    segmentPointStart[0] = grid->nodesMap.at(edgeVertices.first)->GetX();
                    segmentPointStart[1] = grid->nodesMap.at(edgeVertices.first)->GetY();
                    segmentPointEnd[0] =grid->nodesMap.at(edgeVertices.second)->GetX();
                    segmentPointEnd[1] = grid->nodesMap.at(edgeVertices.second)->GetY();
                    if (CheckSegmentsCollision(segmentPointStart, segmentPointEnd, midPoint, tryCoors))
                    {
                        flag = false;
                        break;
                    }
                }
                if (flag)
                {
                    firstTriangleCreated = true;
                    points[2] = it->first;
                    coors[4] = tryCoors[4];
                    coors[5] = tryCoors[5];
                }
            }
        }
    }
    delete[] coors;
    delete[] tryCoors;
    if (firstTriangleCreated)
    {
        return points;
    }
    return nullptr;
}

bool DelaunayTriangulation2D::CheckSegmentsCollision(double*segmentPointStart, double* segmentPointEnd, double* linePointStart, double* tryCoors)
{
    double* linePointEnd{ new double[2] };
    linePointEnd[0] = tryCoors[4];
    linePointEnd[1] = tryCoors[5];

    double intersectPoint[2]{ 0,0 };

    bool isAFunction1{ (abs(segmentPointStart[0] - segmentPointEnd[0]) >= pow(10, -10)) };
    bool isAFunction2{ (abs(linePointStart[0] - linePointEnd[0]) >= pow(10, -10)) };

    if (isAFunction1 && isAFunction2)
    {
        double alpha1{ (segmentPointStart[1] - segmentPointEnd[1]) / (segmentPointStart[0] - segmentPointEnd[0]) };
        double beta1{ segmentPointStart[1] - alpha1 * segmentPointStart[0] };
        double alpha2{ ((linePointEnd[1] - linePointStart[1]) / (linePointEnd[0] - linePointStart[0])) };
        double beta2{ linePointStart[1] - alpha2 * linePointStart[0] };

        intersectPoint[0] = (beta2 - beta1) / (alpha1 - alpha2);
        intersectPoint[1] = alpha1 * intersectPoint[0] + beta1;

    }
    else if(isAFunction1)
    {
        double alpha1{ (segmentPointStart[1] - segmentPointEnd[1]) / (segmentPointStart[0] - segmentPointEnd[0]) };
        double beta1{ segmentPointStart[1] - alpha1 * segmentPointStart[0] };

        intersectPoint[0] = linePointStart[0];
        intersectPoint[1] = alpha1 * intersectPoint[0] + beta1;


    }
    else if (isAFunction2)
    {
        double alpha2{ ((linePointEnd[1] - linePointStart[1]) / (linePointEnd[0] - linePointStart[0])) };
        double beta2{ linePointStart[1] - alpha2 * linePointStart[0] };

        intersectPoint[0] = segmentPointStart[0];
        intersectPoint[1] = alpha2 * intersectPoint[0] + beta2;
    }

    if (intersectPoint[0] >= std::min(segmentPointStart[0], segmentPointEnd[0]) && intersectPoint[0] >= std::min(linePointStart[0], linePointEnd[0])
        && intersectPoint[0] <= std::max(segmentPointStart[0], segmentPointEnd[0]) && intersectPoint[0] <= std::max(linePointStart[0], linePointEnd[0])
        && intersectPoint[1] >= std::min(segmentPointStart[1], segmentPointEnd[1]) && intersectPoint[1] >= std::min(linePointStart[1], linePointEnd[1])
        && intersectPoint[1] <= std::max(segmentPointStart[1], segmentPointEnd[1]) && intersectPoint[1] <= std::max(linePointStart[1], linePointEnd[1]))
    {
        delete[] linePointEnd;
        return true;
    }

    delete[] linePointEnd;
    return false;
}

void DelaunayTriangulation2D::AddTriangle(int* elementID, int* points)
{
    int point0{ points[0] };
    int point1{ points[1] };
    int point2{ points[2] };

    int* pointsIDs{ nullptr };
    pointsIDs = new int[3];
    pointsIDs[0] = point0;
    pointsIDs[1] = point1;
    pointsIDs[2] = point2;
    this->grid->elementsMap[*elementID] = new CST(*elementID, pointsIDs);
    *elementID = *elementID + 1;
}

void DelaunayTriangulation2D::PassSegments(std::vector<std::pair<int, int>>* outerBoundaryEdges, std::vector<std::pair<int, int>>* innerBoundaryEdges)
{
    std::pair<int, int> newEdge{};
    for (const auto edge : *outerBoundaryEdges)
    {
        newEdge = std::make_pair(edge.first, edge.second);
        grid->segments.push_back(newEdge);
    }
    for (const auto edge : *innerBoundaryEdges)
    {
        newEdge = std::make_pair(edge.first, edge.second);
        grid->segments.push_back(newEdge);
    }
}

void DelaunayTriangulation2D::PassVertices(std::map<int, std::pair<double, double>>* vertices)
{
    for (auto &vertice : *vertices)
    {
        grid->nodesMap[vertice.first] = new Node(vertice.first, vertice.second.first, vertice.second.second);
    }
}

double DelaunayTriangulation2D::Orient(double* coors)
{
    return (coors[0] - coors[4]) * (coors[3] - coors[5]) - (coors[1] - coors[5]) * (coors[2] - coors[4]);
}

bool DelaunayTriangulation2D::InsideArea(double* coors)
{
    double value{ 0 };
    double table00{ coors[0] - coors[6]};
    double table01{ coors[1] - coors[7] };
    double table02{ pow(table00, 2) + pow(table01, 2) };
    double table10{ coors[2] - coors[6]};
    double table11{ coors[3] - coors[7]};
    double table12{ pow(table10,2) + pow(table11,2) };
    double table20{ coors[4] - coors[6]};
    double table21{ coors[5] - coors[7]};
    double table22{ pow(table20,2) + pow(table21,2) };

    value = table00 * (table11 * table22 - table12 * table21) - table01 * (table10 * table22 - table12 * table20) + table02 * (table10 * table21 - table11 * table20);
    if (value > -pow(10, -10))
    {
        return true;
    }
    return false;
}