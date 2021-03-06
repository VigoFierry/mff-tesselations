// Delaunay.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <algorithm>  // for min and max
#include <vector>
#include <iterator>
#include <random>
#include <atomic>
#include<unordered_map>


// NOTES
// o Classes in separate files
// o Understand power diagrams and incircle in them
// o Extend the code for power diagrams

// o Walking algorithm
// o Tetrahedron representation


// Abstract class shape that will contain all possible facets
// Useful for e.g. list of a certain k-facet?
class Shape {
public:
	virtual ~Shape() {};
};



class Point{
private:
	double x, y, z;
	int id;

	static std::random_device seed_generator;
	static unsigned seed;
	static std::mt19937 mersenne_generator;
	static std::uniform_real_distribution<double> distribution;

public:	
	static std::atomic<int> s_id;

	Point(double _x = 0.0, double _y = 0.0, double _z = 0.0) : x(_x), y(_y), z(_z), id(s_id) { 
		s_id++; 
	};
	
	// Point(const Point& other) : x(other.getX()), y(other.getY()), z(other.getZ) {};   No need for a copy constructor - it is constructed automatically
	
	static Point randomizePoint() {
		return Point(distribution(mersenne_generator), distribution(mersenne_generator), distribution(mersenne_generator));
	};



	void setX(double val) { x = val; };
	void setY(double val) { y = val; };
	void setZ(double val) { z = val; };

	double getX() const { return x; };
	double getY() const { return y; };
	double getZ() const { return z; };
	int getid() const { return id; };

	Point operator+(const Point& p) {
		return Point(x + p.x, y + p.y, z + p.z);
	};
	Point operator-() {
		return Point(-x, -y, -z);
	};

	Point operator*(double factor) {
		return Point(x*factor, y*factor, z*factor);
	};

	double operator*(Point p) {
		return x*p.getX() + y*p.getY() + z*p.getZ();
	};

	Point operator=(const Point& rhs) {
		if (this == &rhs)
			return *this;
		this->setX(rhs.getX());
		this->setY(rhs.getY());
		this->setZ(rhs.getZ());
		return *this;
	};

	double Length() {
		return sqrt((*this)*(*this));
	};

	bool operator==(const Point& rhs) {
		return (this->getX() == rhs.getX() && this->getY() == rhs.getY() && this->getZ() == rhs.getZ());
	};

};

// Random
std::random_device Point::seed_generator;
unsigned Point::seed = seed_generator();
std::uniform_real_distribution<double> Point::distribution(0, 1);
std::mt19937 Point::mersenne_generator(Point::seed);

// Initialization of the atomic static counter - counts points
std::atomic<int> Point::s_id = 0;


Point operator-(Point p, Point q) {
	return p + (-q);
};

Point operator*(double factor, Point p) {
	return p*factor;
};


std::ostream& operator<<(std::ostream& out, Point p) {
	out << "(" << p.getX() << "," << p.getY() << "," << p.getZ() << ")";
	return out;
};

class Line{
private:
	Point a, b;

public:
	Line(Point _a, Point _b) : a(_a), b(_b) {};

	void setA(Point p) { a = p; };
	void setB(Point p) { b = p; };

	Point getA() const { return a; };
	Point getB() const { return b; };

	double Length() const { return sqrt((b - a)*(b - a)); };
};


std::ostream& operator<<(std::ostream& out, Line l) {
	out << "[" << l.getA() << "," << l.getB() << "]";
	return out;
};

class Face{
private:
	Point a, b, c;
public:
	Face(Point _a, Point _b, Point _c) : a(_a), b(_b), c(_c) {};
	Face(std::vector<Point> _v) : a(_v[0]), b(_v[1]), c(_v[2]) {};
	Face() : a(), b(), c() {}; // Default constructor for default HashNode

	void setA(Point p) { a = p; };
	void setB(Point p) { b = p; };
	void setC(Point p) { c = p; };

	Point getA() const { return a; };
	Point getB() const { return b; };
	Point getC() const { return c; };

	std::vector<Point> getPoints() const { 
		std::vector<Point> v = { a,b,c }; 
		return v;
	};

	
	double minEdgeLength() {
		return std::min((b - a).Length(), std::min((c - a).Length(), (c - b).Length()));
	};

	// https://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
	double area() {
		Point AB = b - a;
		Point AC = c - a;
		
		return (1 / 2.0 * sqrt(pow(AB.getY()*AC.getZ() - AB.getZ()*AC.getY(), 2) + pow(AB.getZ()*AC.getX()-AB.getX()*AC.getZ(),2) + pow(AB.getX()*AC.getY() - AB.getY()*AC.getX(),2) ) );
	};


	bool operator==(const Face& rhs) {
		return (this->getA() == rhs.getA() && this->getB() == rhs.getB() && this->getC() == rhs.getC());
	};
};



bool operator==(const Face& lhs, Face rhs) {
	return rhs == lhs;
};

bool operator!=(Face f, Face g) {
	return !(f == g);
};

std::ostream& operator<<(std::ostream& out, Face f) {
	out << "[" << f.getA() << "," << f.getB() << f.getC() << "]";
	return out;
};



class Tetrahedron {
private:
	Point a, b, c, d;
public:
	Tetrahedron(Point _a, Point _b, Point _c, Point _d) : a(_a), b(_b), c(_c), d(_d) {};
	Tetrahedron(Face _f, Point _d) : a(_f.getA()), b(_f.getB()), c(_f.getC()), d(_d) {};
	Tetrahedron(std::vector<Point> _v) : a(_v[0]), b(_v[1]), c(_v[2]), d(_v[3]) {};

	void setA(Point p) { a = p; };
	void setB(Point p) { b = p; };
	void setC(Point p) { c = p; };
	void setD(Point p) { d = p; };

	Point getA() const { return a; };
	Point getB() const { return b; };
	Point getC() const { return c; };
	Point getD() const { return d; };

};





void getCofactor(double m[5][5], double cofactor[5][5], int p, int q, int n) {
	int i = 0, j = 0;

	// Loop over matrix m
	for (int row = 0; row < n; row++) {
		for (int col = 0; col < n; col++) {

			if (row != p && col != q) {
				cofactor[i][j++] = m[row][col];
			};

		};
		// end of row
		if (j != 0) { j = 0; i++; };
	};
};

void writeMatrix(double m[5][5], int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << m[i][j] << " ";
		};
		std::cout << std::endl;
	};
};

double det(double m[5][5], int n) {
	double D = 0.0;

	if (n == 2) {
		return m[0][0] * m[1][1] - m[1][0] * m[0][1];
	}
	else {
		double cofactor[5][5];
		int sign = 1;
		// Iterate over the first row
		for (int i = 0; i < n;i++) {
			// Get cofactor matrix

			getCofactor(m, cofactor, 0, i, n);
			D += sign*m[0][i] * det(cofactor, n - 1);
			sign = -sign; // Sign changes each time
		};
		return D;
	};
};


void INCIRCLE(Point p0, Point p1, Point p2, Point p3, Point p4) {
	double m[5][5];

	// Weights. For power diagram would also include the radii
	double w0 = pow(p0.getX(), 2) + pow(p0.getY(), 2) + pow(p0.getZ(), 2);
	double w1 = pow(p1.getX(), 2) + pow(p1.getY(), 2) + pow(p1.getZ(), 2);
	double w2 = pow(p2.getX(), 2) + pow(p2.getY(), 2) + pow(p2.getZ(), 2);
	double w3 = pow(p3.getX(), 2) + pow(p3.getY(), 2) + pow(p3.getZ(), 2);
	double w4 = pow(p4.getX(), 2) + pow(p4.getY(), 2) + pow(p4.getZ(), 2);

	// Define INCIRCLE matrix
	m[0][0] = p0.getX();  m[0][1] = p0.getY(); m[0][2] = p0.getZ(); m[0][3] = w0; m[0][4] = 1;
	m[1][0] = p1.getX();  m[1][1] = p1.getY(); m[1][2] = p1.getZ(); m[1][3] = w1; m[1][4] = 1;
	m[2][0] = p2.getX();  m[2][1] = p2.getY(); m[2][2] = p2.getZ(); m[2][3] = w2; m[2][4] = 1;
	m[3][0] = p3.getX();  m[3][1] = p3.getY(); m[3][2] = p3.getZ(); m[3][3] = w3; m[3][4] = 1;
	m[4][0] = p4.getX();  m[4][1] = p4.getY(); m[4][2] = p4.getZ(); m[4][3] = w4; m[4][4] = 1;


	// Calculate the determinant of the INCRICLE matrix
	std::cout << "Det " << det(m, 5) << std::endl;

	// writeMatrix(m, 5);	
};

void copyMatrix(double A[5][5], double B[5][5], int n) {
	for (int i = 0;i < n;i++) {
		for (int j = 0; j < n; j++) {
			B[i][j] = A[i][j];
		};
	};
};




std::vector<double> getBarycentricCoordinates(Point p0, Point p1, Point p2, Point p3, Point p) {
	double m[5][5];

	m[0][0] = p0.getX(); m[0][1] = p1.getX(); m[0][2] = p2.getX(); m[0][3] = p3.getX();
	m[1][0] = p0.getY(); m[1][1] = p1.getY(); m[1][2] = p2.getY(); m[1][3] = p3.getY();
	m[2][0] = p0.getZ(); m[2][1] = p1.getZ(); m[2][2] = p2.getZ(); m[2][3] = p3.getZ();
	m[3][0] = 1;		 m[3][1] = 1;		  m[3][2] = 1;		   m[3][3] = 1;

	// std::cout << det(m, 4) << std::endl;

	std::vector<double> coordinates;
	double mCram[5][5];
	double det_m = det(m, 4);
	

	for (int i = 0; i < 4; i++) {
		copyMatrix(m, mCram, 4);
		mCram[0][i] = p.getX();
		mCram[1][i] = p.getY();
		mCram[2][i] = p.getZ();
		mCram[3][i] = 1;

		coordinates.push_back(det(mCram,4)/det_m);
	}


	return coordinates;
};

std::vector<double> getBarycentricCoordinates(Tetrahedron tetra, Point p) {
	return getBarycentricCoordinates(tetra.getA(), tetra.getB(), tetra.getC(), tetra.getD(), p);
};

// Computes barycentric coordinates of point p with respect to the tetrahedron formed by face and v
std::vector<double> getBarycentricCoordinates(Face face, Point v, Point p) {
	return getBarycentricCoordinates(face.getA(), face.getB(), face.getC(), v, p);
};





// Learn about: templates, iterators, copy
template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
	if (!v.empty()) {
		out << '[';
		std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
		out << "\b\b]";
	}
	return out;
}




int pokusna(Point& p) {
	return p.getid();
}

// Finds the index of element e in vector v
// Why not const std::vector<T>& v
template <typename T>
int vectorIndex(const T& e, std::vector<T>& v) {
	return std::distance(v.begin(), (std::find(v.begin(), v.end(), e)));
}




int testhash(const Face& face) {
	const unsigned int M = 7;
	int n = 100;
	int a = (face.getA()).getid();
	int b = (face.getB()).getid();
	int c = (face.getC()).getid();
	return (a*n*n + b * n + c) % M;
}

// Hash function for unordered map for Face

struct faceHash {
	size_t operator()(const Face& face) const {
		const unsigned int M = 10007;
		int n = 100;
		int a = (face.getA()).getid();
		int b = (face.getB()).getid();
		int c = (face.getC()).getid();
		return (a*n*n + b * n + c) % M;
	}
};


// typedef std::pair <int, int> intpair;
class HashNode {
private:
	Face face;
	int e;
	int f;
	// HashNode* next; Not needed since unordered map resolves colisions by itself?

public:
	HashNode(const Face &_face, int _e = -1, int _f = -1) : face(_face), e(_e), f(_f) {};
	HashNode() : face(), e(-1), f(-1) {};


	Face getFace() const {
		return face;
	};

	void setFace(Face face) {
		HashNode::face = face;
	};

	int getPointe() const {
		return HashNode::e;
	};

	int getPointf() const {
		return HashNode::f;
	};

	void setPointe(int e) {
		HashNode::e = e;
	};

	void setPointf(int f) {
		HashNode::f = f;
	};

	// HashNode* getNext() const {
	//	return next;
	// };


	// void setNext(HashNode* next) {
	//	HashNode::next = next;
	// };

};


std::ostream& operator<<(std::ostream& out, HashNode h) {
	out << h.getFace() << "," << h.getPointe() << "," << h.getPointf();
	return out;
};






void Incremental() {
	std::vector<Point> Vertices;

	Point p(1,1,2);
	Point q(0, 0, 1);
	Vertices.push_back(q);
	Vertices.push_back(p);
	std::cout << Vertices << Vertices[0] << std::endl;
	std::cout << vectorIndex(p, Vertices) << std::endl;	
};

// Template od Nohice
template <typename T>
using matrix = std::vector<std::vector<T>>;


/*	Walking algorithm
	int passingFace() { returns 0 if inside, 1-4 for number of a face through which to walk - but stochastic?}
		Where to walk? Minus sign means walk on the side where there's a zero

	passThroughFace() { Returns tetrahedron (?)  }
	main() {
		Pick a random tetrahedron (how?)
		passingFace()
	}
	


*/

int main()
{
		
	
	
	
	// Initialize a point
	Point point(1, 2, 3);
	// Initialize default point
	Point zero;
	// Possible operations with a point
	// Each such operation initializes a new point with the value of the expression
	std::cout << point << std::endl;
	point + point;
	(-point);
	2 * point;
	point * 2;
	point.setX(10);
	
	// Initializing a line;
	Line line(point, zero);
	std::cout << line << std::endl;
	std::cout << line.Length() << std::endl;

	// Testing INCIRCLE
	// Define points for testing incircle
	Point p0(1, 0, 0);
	Point p1(-1, 0, 0);
	Point p2(0, 1, 0);
	Point p3(0, -1, 0);
	Point p4(0, 0, 0.9);

	INCIRCLE(p0, p1, p2, p3, p4);


	// Barycentric cooods of multiple points
	//std::vector<Point> Vertices;

	//for (int i = 0; i <= 100; i++) {
	//	Vertices.push_back(Point(true));
	//};

	//std::cout << Vertices << std::endl;

	//matrix<double> MapVertices;

	//std::transform(Vertices.begin(), Vertices.end(), std::back_inserter(MapVertices), [&](Point p) { return getBarycentricCoordinates(r0, r1, r2, r3, p); });
	////std::cout << MapVertices << std::endl;

	//for (int i = 0; i < MapVertices.size(); i++) {
	//	std::cout << Vertices[i] << ' ' << MapVertices[i] << std::endl;
	//};
	//
	
	// Nohic solution to a random point
	Point pokusnyBod = Point::randomizePoint();
	std::cout << pokusnyBod << std::endl;




	//// Testing basics of incremental
	// Create some points
	std::vector<Point> Vertices;
	// Main face
	Vertices.push_back(Point(0, 0, 0));
	Vertices.push_back(Point(1, 0, 0));
	Vertices.push_back(Point(0.5, 1, 0));
	// And two points creating two tetrahedrons
	Vertices.push_back(Point(0.5, 0.5, 1));
	Vertices.push_back(Point(0.5, 0.5, -1));
	// Check if all ok
	std::cout << Vertices << std::endl;
	std::cout << Vertices[4] << std::endl;
	// Create the basic face and its hashnode
	Face face(Vertices[0], Vertices[1], Vertices[2]);
	HashNode hnode(face);
	hnode.setPointe(3);
	hnode.setPointf(4);
	// Create HashMap
	std::unordered_map<Face, HashNode, faceHash> faces;
	faces.insert(std::make_pair(face, hnode));
	std::cout << faces[face] << std::endl;	
	// Create a point that's inside the upper tetrahedron
	Vertices.push_back(Point(0.5, 0.5, 0.1));
	// Start in lower tetrahedron and check where to go
	std::vector<Point> tetra(face.getPoints());
	tetra.push_back(Vertices[4]);
	std::vector<double> result = getBarycentricCoordinates(Tetrahedron(tetra), Vertices[5]);
	int sideIndex = std::distance(result.begin(),std::find_if(result.begin(), result.end(), [](int i) { return i < 0; })) - 1;
	std::cout << tetra[sideIndex] << std::endl;
		

	

	
	// New HashMap === std::unordered_map
	
	














	system("pause");

	return 0;
}

