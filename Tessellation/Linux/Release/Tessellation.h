#include <vector>
#include <queue>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <math.h>
#include <memory>
#include <algorithm>
#include <string>
#include <cstring>
using namespace std;

/* Author: Milan Pultar, milan.pultar@gmail.com, date: 2017-11-06*/
/* One can modify classes below by inheritance (feel free to use virtual methods where needed); you can't modify the functions itself */
/* BEWARE: when writing new functions, always keep in mind that functions with raw pointer as an argument should return quickly, they should not call other methods.
	conversion shared_ptr<> -> *[raw] is considered unsafe */

#define pow2(x) (x)*(x)
#define pow3(x) (x)*(x)*(x)
// ^ macros

class Tessellation;
class Vertex;
class Edge;
class Face;
class Grain;
class Plane;
// ^ forward declaration

class Tessellation {
		// ATTRIBUTES:
public:	vector<shared_ptr<Grain>> Grains;
		// ^ these grains are in the tesselation, use this to access the geometry
		inline size_t get_number_of_grains() { return Grains.size();  };
		// get the number of grains in this tessellation
		vector<shared_ptr<Face>> WindowFaces;
		// ^ These are the faces of the current window. For example, this window could be created by the default constructor.
		// METHODS:
		vector<shared_ptr<Face>> GetWindow();
		// ^ Creates a deep copy (!) of WindowFaces and returns it.
		void UpdateWindowFaces(Plane*, vector<shared_ptr<Face>>&, Grain*, Grain*);
		// ^ update window faces, which are in the vector<shared_ptr<Face>>&, by the plane
		shared_ptr<Grain> GetClosestGrain(Vertex*);
		// ^ gets the closest grain in this->Grains to the Vertex*, measured with Grain::DistanceTo method
		void PrintEdges(string, bool only_once = false);
		// ^ Prints the tesselation to output.
		// ^ The format is: On each line there are 6 numbers: (x y z) (x y z) describing a line. Each grain structure is ended by  NaN 6 times. 
		// ^ At the end after 2 lines o NaN the centers of grains are listed. The format is (x y z).
		// ^ the bool only_once indicates whether we want to print each vertex only once when the Grains are in the connected form (e.g. Grains.back() = Grains.back().to_connected())
		int InsertGrain(double, double, double, double);
		// ^ Use this method to insert standard Grain, arguments: x,y,z,radius; returns false if singularity detected. This is very rare, it prevents the tessellation from corruption.
		// ^ return value: true if inserted; false if corruption was detected, insertion procedure was blocked.
		int InsertGrain(shared_ptr<Grain>);
		// ^ Use this method to insert Grain. You may insert your own modified grain as well (It should be a child of the standard Grain of course).
		bool RemoveGrain(shared_ptr<Grain>);
		// ^ Use this method to remove a grain.
		// ^ Argument: the pointer to the grain to remove. It is your responsibility to pass a correct argument (e.g. Tessellation->Grains.back()) ! Method will fail otherwise.
		// ^ return value: true if removed; false if corruption was detected, remove procedure was blocked. Happens very rarely.
		// CONSTRUCTORS:
		Tessellation();
		// ^ Creates default window, it has coord. of vertices in {0,1}
		// ^ You can create another constructor to suit your needs (or you can just modify the faces before adding any grains).
		// ^ Tessellation class should do the work right for any convex window.
};
class Vertex {
		// ATTRIBUTES:
public: double x, y, z;
		// METHODS
		Grain* WhoIsCloser(Grain*, Grain*);
		// ^ Returns a pointer to the grain which is closer to (Vertex)this; measured by (*Grain).DistanceTo() method.
		Vertex CrossProduct(Vertex*);
		// ^ Returns a cross product of (Vertex)this and the Vertex in the argument (Vertex is considered as a vector here).
		size_t GetClosest(vector<shared_ptr<Vertex>>);
		// ^ Returns the index of the closest Vertex in vector in the argument to (Vertex)this. It skips all pointers equal to (Vertex)this
		double EuclidDist(Vertex*);
		// ^ Computes Euclidean distance between Vertex in the argument and (Vertex)this.
		// ^ Returns the Euclid. distance.
		// CONS, DESC, OPERATORS:
		Vertex(double Ix, double Iy, double Iz) : x(Ix), y(Iy), z(Iz) {};
		//Vertex::Vertex(Grain*);
		Vertex(Vertex* IB) : x(IB->x), y(IB->y), z(IB->z) {};
		// ^ simple constructors for ease of use
		virtual Vertex operator/(const double num) { return Vertex(this->x / num, this->y / num, this->z / num); };
		virtual inline bool operator==(const Vertex &IA) { return IA.x == this->x && IA.y == this->y && IA.z == this->z; }
		virtual inline Vertex operator-() { return Vertex(-this->x, -this->y, -this->z); }
		virtual inline Vertex operator+(const Vertex &IB) { return Vertex(this->x + IB.x, this->y + IB.y, this->z + IB.z); }
		// ^ defines scalar*vector operations and unary operations
};
class Edge {
		// ATTRIBUTES:
public:	shared_ptr<Vertex> u, v;
		// Edge(this) connects the (Vertex)u and (Vertex)v.
		// METHODS:
		bool IsOnOneSide(Grain* A, Grain* B);
		// ^ Returns bool indicating whether both (Vertex)u and (Vertex)v are closer to (Grain)A.
		// ^ Use this method to indicate whether the edge has been destroyed while inserting new Grain (only if you intend to modify the insertion procedure).
		inline double EuclidDist(Edge*);
		// ^ Computes Euclidean distance between Edge in the argument and (Edge)this. Returns the max distance between all 4 pairs of vertices.
		// CONSTRUCTORS, DESCTRUCTORS
		Edge(shared_ptr<Vertex> Iu, shared_ptr<Vertex> Iv) : u(Iu), v(Iv) {};
		Edge() {};
		// ^ (Vertex)u and (Vertex)v have to be set before using any method in Edge class !
};
class Face {
		// ATTRIBUTES:
public: shared_ptr<Grain> a, b;
		// ^ the two Grains between which (Face)this lies, b is NULL if this is a part of the window (no outer grain)
		shared_ptr<Plane> plane = NULL;
		// ^ what plane is defined by this face, not always evaluated, therefore may be NULL
		// METHODS:
		vector<shared_ptr<Edge>> Edges;
		// ^ Edges of this
		bool CutByPlane(Plane*, Grain*, Grain*);
		// ^ cuts (Face)this by (Plane)pl given the two Grains in,out. It must hold that (Plane)pl is the set of all points which have the same distance to (Grain)in, (Grain)out
		// ^ returns false if (Face)this lies completely in (Grain)out
		bool lies_in_first(Grain*, Grain*);
		// ^ whether all vertices of the Face lie closer to the FIRST Grain
		virtual void set_plane(Plane *Ir) { this->plane = shared_ptr<Plane>(Ir); };
		// CONS, DESC:
		Face() {};
		// ^ (Face)a and (Face)b have to be set before using any method in Face class !
		Face(Plane *Ir) { this->plane = shared_ptr<Plane>(Ir); };
		// ^ This is a dangerous constructor, use only if you understand basics of shared_ptr<>. The argument Ir is a unique raw pointer to Plane instance, which was probably just created.
};
class Grain {
		// ATTRIBUTES
public:	double x, y, z, r;
		// [x,y,z] is the center of the Grain. r is the radius of the sphere of the grain.
		vector<shared_ptr<Face>> Faces;
		// ^ faces of the grain, unconnected representation
		vector<shared_ptr<Grain>> AdjacentGrains;
		// ^ Pointers to all adjacent grains.
		inline vector<shared_ptr<Grain>> get_adjacent_grains() { return AdjacentGrains;  };
		// ^ Creates a deep copy of the vector of AdjacentGrains. Use this to get all the neighbouring Grains.
		// METHODS
		virtual double DistanceTo(Vertex*);
		// ^ uses power distance ; other distances work as well if this method and method Plane->ComputeD are changed accordingly (under certain presumptions...)
		double VolumeApprox(double precision = 0.05);
		// ^ DEPRECATED
		// ^ computes volume, but is very inaccurate, analytic solution recommended.
		bool OnSameSide(Vertex*, Vertex*);
		// ^ Returns bool indicating whether both points in arguments lie in the (Grain)this.
		shared_ptr<Grain> to_connected();
		// ^ Returns Grain, where each object exists only once
		Vertex get_center();
		// ^ returns Vertex(x,y,z)
		// CONSTRUCTORS, DESCTRUCTORS
		Grain(double Ix, double Iy, double Iz, double Ir) : x(Ix), y(Iy), z(Iz), r(Ir) { };
		// ^ creates a grain with center Ix,Iy,Iz, radius Ir
		Grain() { };
		// ^ creates a grain
};
class Plane {
		// ATTRIBUTES
public:	double a, b, c, d;
		// ^ format: a*x + b*y + c*z = d
		// METHODS
		virtual char WhichSide(double Ix, double Iy, double Iz) { return (a*Ix + b*Iy + c*Iz > d) ? 1 : 0; };
		virtual char WhichSide(Vertex *Ib) { return this->WhichSide(Ib->x, Ib->y, Ib->z); };
		virtual bool IntersectsEdge(shared_ptr<Edge> In) { return WhichSide(In->u->x, In->u->y, In->u->z) != WhichSide(In->v->x, In->v->y, In->v->z); };
		// ^ Returns bool indicating whether (Plane)this intersects (Edge)In in the argument.
		shared_ptr<Vertex> IntersectWithLine(Vertex*, Vertex*);
		// ^ Returns the intersection Vertex where the line connecting two Vertices in the argument meets the (Plane)this.
		void ComputeD(Grain*, Grain*);
		// ^ Computes the number d, given the two grains. This method is written specially for the power distance.
		// ^ Override this method if you wish to create a tessellation based on another distance function !
		virtual void ComputeD(Vertex* IB) { this->d = this->a * IB->x + this->b * IB->y + this->c * IB->z; };
		// ^ Computes the d, given (Vertex)IB lying in the plane.
		double DistanceTo(Vertex*);
		// ^ Computes the smallest distance (Euclidean) from the Vertex in the argument to the (Plane)this
		// CONS, DESC:
		Plane(Grain*, Grain*);
		// ^ computes a,b,c, given two grains, then calls ComputeD with the same args (so you don't have to call it yourself)
		Plane(double Ia, double Ib, double Ic) : a(Ia), b(Ib), c(Ic) {};
		// ^ construct plane if a,b,c, are known
};
