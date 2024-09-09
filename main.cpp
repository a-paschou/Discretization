//#include <vector>
//#include <iostream>
//#include "Surface.h"
//#include "Line.h"
//#include "Circle.h"
//#include "PoissonDiskSampling.h"
//#include "DelaunayTriangulation2D.h"
//#include <iostream>
//#include <fstream>
//
//void WriteMesh(Grid*);
//void main()
//{
//	IShape* A = new Point(0, 0);
//	IShape* B = new Point(170, 0);
//	IShape* C = new Point(170, 100);
//	IShape* D = new Point(0, 100);
//
//	IShape* AB = new Line((Point*)A, (Point*)B);
//	IShape* BC = new Line((Point*)B, (Point*)C);
//	IShape* CD = new Line((Point*)C, (Point*)D);
//	IShape* DA = new Line((Point*)D, (Point*)A);
//
//	Surface surface{ std::vector<IShape*> { AB, BC, CD, DA } };
//
//	IShape* circle1 = new Circle(new Point(85, 50), 25);
//	//IShape* circle2 = new Circle(new Point(3, 3), 2);
//
//	
//	
//	//IShape* cutPoint1 = new Point(20, 20);
//	//IShape* cutPoint2 = new Point(40, 20);
//	//IShape* cutPoint3 = new Point(40, 40);
//	//IShape* cutPoint4 = new Point(20, 40);
//
//	//IShape* cutEdge1 = new Line((Point*)cutPoint1, (Point*)cutPoint2);
//	//IShape* cutEdge2 = new Line((Point*)cutPoint2, (Point*)cutPoint3);
//	//IShape* cutEdge3 = new Line((Point*)cutPoint3, (Point*)cutPoint4);
//	//IShape* cutEdge4 = new Line((Point*)cutPoint4, (Point*)cutPoint1);
//	//std::vector<IShape*> cutQuad{ cutEdge1, cutEdge2, cutEdge3, cutEdge4 };
//
//	
//	surface.Cut(std::vector<IShape*> {circle1});
//	//surface.Cut(cutQuad);
//
//	double radius{ 2.5 };
//	int k{ 20 };
//
//	PoissonDiskSampling test{ &surface, radius, k };
//	IDelaunayTriangulation* mesh = new DelaunayTriangulation2D(&test.vertices, &test.outerBoundaryEdges2, &test.innerBoundaryEdges2);
//	
//	WriteMesh(&mesh->grid);
//
//	std::cout << "End" << std::endl;
//}
//
//void WriteMesh(Grid* grid)
//{
//	int elementType{ 2 };
//	std::ofstream MyFile("mesh.txt");
//
//	MyFile << "$MeshFormat";
//	MyFile << "\n4.1 0 8 MSH4.1, ASCII";
//	MyFile << "\n$EndMeshFormat";
//
//	MyFile << "\n$Nodes";
//	MyFile << "\n1 " << grid->vertices->size() << " 1 " << grid->vertices->size();
//	MyFile << "\n2 1 0 " << grid->vertices->size();
//
//	for (std::map<int, std::pair<double, double>>::iterator it = grid->vertices->begin(); it != grid->vertices->end(); ++it)
//	{
//		MyFile << "\n" << (it->first + 1);
//	}
//	for (std::map<int, std::pair<double, double>>::iterator it = grid->vertices->begin(); it != grid->vertices->end(); ++it)
//	{
//		MyFile << "\n" << it->second.first << " " << it->second.second << " 0";
//	}
//	MyFile << "\n$EndNodes";
//
//	MyFile << "\n$Elements";
//	MyFile << "\n1 " << grid->elements.size() << " 1 " << grid->elements.size();
//	MyFile << "\n2 1 " << elementType << " " << grid->elements.size();
//	for (std::map<int, int*>::iterator it = grid->elements.begin(); it != grid->elements.end(); it++)
//	{
//		MyFile << "\n" << (it->first + 1) << " " << (it->second[0] + 1) << " " << (it->second[1] + 1) << " " << (it->second[2] + 1);
//	}
//	MyFile << "\n$EndElements";
//
//	MyFile.close();
//}