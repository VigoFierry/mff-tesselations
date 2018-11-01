#pragma once

#include <iostream>
#include <fstream>
#include <random>  // package for generating random variables - http://en.cppreference.com/w/cpp/numeric/random or http://itscompiling.eu/2016/04/11/generating-random-numbers-cpp/
#include <vector>  // vectors

#include <cstdio>  // work with files
#include <math.h>  // math functions: tan, pow, sqrt, ...
#include <time.h>  // time 

#include <omp.h>   // OpenMP for parallel programming

// Voro++ should be included in directories
#include "voro++.hh"  // Voro++ - http://math.lbl.gov/voro++/

// Eigen should be included in directories
//#include <Eigen/Core> // dependency of LBFGS++, http://eigen.tuxfamily.org/index.php?title=Main_Page

// LBFGS++ should be included in directories 
//#include <LBFGS.h>  // LBFGS optimalization algorithm, https://yixuan.cos.name/LBFGSpp/


/* Author: Filip Seitl, seitlf@seznam.cz, date: 2017-12-20 */
/* uncompleted version */

/* This version is not final, a lot of functions is under development. To avoid undesirable errors is recomended to follow the example and use mostly
   the functions mentioned there. Unstable can be parts concerning Laguerre tessellation and different types of energy function than is V2, as well the end
   of estimation section was not tested at all. 
   
   A lot of notes in cpp files is for need of author and is irrelevant for the final implementation. */




// prob distributions //
// - fcs generating random variables
double rnd();
double uniform(int, int);
int uniform_int(int a, int b);
double normal(double mu, double sigma);
double poisson(int);
double gamma(double alfa, double beta);
double triangle(double a, double b, double c);
double exponential(double lambda);



// numeric solution of equation //
double f(double x, std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int &k, int &y);
// ^ -- definition of the function whose root should be found; in the pseudolikelihood case it is derivative of PLL function with respect to theta,
// where intensity parameter zet is expressed from the derivative of PLL w.r.t. zet (PLL = pseudo-loglikelihood contrast function);
// input parameters are obtained from function estim (section estimation below)
double fprime(double x, std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int &k, int &y);
// ^ -- the derivative of function f above used in Newton-Raphson algorithm
double bisection(std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int &k, int y);
// ^ -- bisection method for computing roots of function f, returns the root
double secant(std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int &k, int y);
// ^ -- secant method 
double NR(std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int &k, int y);
// ^ -- Newton-Raphson method



// data_creator //
// - fcs for creating testing datasets and outputing generators of tessellation
void cube_net(double h);
// ^ -- creates a list of generators (form: id x-coordinate y-coordinate z-coordinate) of cubic lattice with cubic edge of length h and starting point [h/2,h/2,h/2] (suitable for Voronoi)
void cube_rad_net(double h);
// ^ -- creates a list where are added radii (form: id x-coordinate y-coordinate z-coordinate radius), all radii are equal to one (suitable for Laguerre)
void random_net(double h);
void random_rad_net(double h);
void random_container(int n);
void write_boxes(bool soubor, voro::container &con);
void write_boxes(bool soubor, voro::container_poly &con);
void write_container(voro::container &con);
// ^ -- the container is written into txt file in the form id x-coordinate y-coordinate z-coordinate; one line = one generator
void write_container(voro::container_poly &con);
// ^ -- the container is written into txt file in the form id x-coordinate y-coordinate z-coordinate radius; one line = one generator
void transform();


// helping_fcs //
// - "small" functions doing smaller calculations for another functions typically
inline int step_int(double a) { return a<0 ? int(a) - 1 : int(a); }

//void positions(std::vector<float> &v);   // nepouzita
void find_pos(int &rijk, int &rq, const int &rid, const voro::container * const pcon);
void find_pos(int &rijk, int &rq, const int &rid, const voro::container_poly * const pcon);
bool find_part(int &ijk, int &q, const int &no, voro::container * const con);
bool find_part(int &ijk, int &q, const int &no, voro::container_poly * const con);

bool erase(const int &rijk, int &rq, std::vector<int> *pfid, voro::container *pcon);
bool erase(const int &rijk, int &rq, std::vector<int> *pfid, voro::container_poly *pcon);
bool erase(const int &rijk, int &rq, voro::container *pcon);
bool erase(const int &rijk, int &rq, voro::container_poly *pcon);
// ^ -- functions for deleting a generator from the container structure (generator is saved in container using two values [ijk,q])

bool are_neighbors(voro::voronoicell_neighbor &rc, const int &rijk, const int &rq, const voro::container * const pcon);
bool are_neighbors(voro::voronoicell_neighbor &rc, const int &rijk, const int &rq, const voro::container_poly * const pcon);
bool secondary(const int cid, const std::vector<int> sec, const voro::container * const pcon, int &ijk, int &q);  // nepouzita
bool terciary(const int cid, const int pid, const std::vector<int> sec, const voro::container * const pcon);
bool terciary(const int cid, const int pid, const std::vector<int> sec, const voro::container_poly * const pcon);
bool terciary(const int cid, const std::vector<int> sec, const voro::container * const pcon);
bool terciary(const int cid, const std::vector<int> sec, const voro::container_poly * const pcon);

bool identical(const std::vector<int> &ra, const std::vector<int> &rb);    // nepouzita
void merge(std::vector<int> &rna, std::vector<int> &rnb, std::vector<int> &ra, std::vector<int> &rb, int k,
	std::vector<int> &raa, std::vector<int> &rbb, std::vector<double> &raaa, std::vector<double> &rbbb);
void merge(std::vector<int> &rna, std::vector<int> &rnb);
bool merge(std::vector<int> &rna, int &rnb);
void vector_dif(std::vector<int> &rna, std::vector<int> &rnb);
void vector_dif(std::vector<int> &rna, std::vector<int> &rnb, std::vector<int> &rnc);

void min_max(double &ra, double &rb);
void min_max(int &ra, int &rb);
void min_max(double &ra, double &rb, double &rc);
void min_max(int &ra, int &rb, int &rc);
// ^ -- returns ordered values

double abs_val(double val);

bool barycentrum(double &rx, double &ry, double &rz, double &rnx, double &rny, double &rnz, double ra, double rb);
// ^ -- returns true iff the barycenter of two cells (given by its centroids/barycenters x,y,z , nx,ny,nz and its volumes) lies within the window 
void bar_coor(double &rx, double &ry, double &rz, double &rnx, double &rny, double &rnz, double &ra, double &rb);
void face_dist(const unsigned int &rfing, const std::vector<int> &fvert, const double &rx, const double &ry, const double &rz,
	double &rnx, double &rny, double &rnz, voro::voronoicell_neighbor &rc);
// ^ -- function computing the normal vector for a given face pointing outside the cell and whose length is equal to the distance of the face to the generator ;
void face_dist(const unsigned int &rfing, const std::vector<int> &fvert, const double &rx, const double &ry, const double &rz,
	const double &tx, const double &ty, const double &tz, double &rnx, double &rny, double &rnz, voro::voronoicell_neighbor &rc);
// ^ -- function computing the normal vector from the barycenter to the given face
void real_coo(double x, double y, double z, double &xx, double &yy, double &zz);
double point_dist(double &rx, double &ry, double &rz, double &rnx, double &rny, double &rnz);
double h_maximum(unsigned int &n, std::vector<int> &vert, voro::voronoicell_neighbor &d, double &x, double &y, double &z);
double h_minimum(unsigned int &n, std::vector<int> &vert, voro::voronoicell_neighbor &d, double &x, double &y, double &z); // nepouzita
void h_fcs(voro::voronoicell_neighbor &d, double &x, double &y, double &z, double &xb, double &yb, double &zb, double &h_max, double &h_min);
void volume_min_max(voro::container &rcon);
void area_min_max(voro::container &rcon);

int find_index(std::vector<int> rv, int &ri);
int common_edge(voro::voronoicell_neighbor &rc, int &f1, int &f2);

void point_density(std::vector<int> counts, voro::container &rcon, int &gsi);
void point_density(std::vector<int> counts, voro::container_poly &rcon, int &gsi);

double ave_rad(voro::container_poly &con);

void new_con(voro::container_poly &con, bool a);

int null_boxes(voro::container &con);
int null_boxes(voro::container_poly &con);

int empty_cells(voro::container_poly &con);
int nonempty_cells(voro::container_poly &con);
void un_vertices(voro::container &rcon);
void test_ed(voro::container &con);

void find05(voro::container_poly &con);




// birth-death-move algorithm (bdma) //
// - fcs controlling the running of BDM-algorithm
bool feasibility(voro::container &con, const double &alfa, const double &beta, const double &B);
bool feasibility(voro::container_poly &con, const double &alfa, const double &beta, const double &B);
// ^ -- returns true if the whole configuration (read tessellation) is feasible w.r.t. hardcore parameters alfa, beta, B
bool overlap_f(voro::container_poly &con, const double &iota);

bool feasibility_subwindows(voro::container_poly &con, const double &alfa, const double &beta, const double &B, double x, double y, double z, double rad);

bool feasibility(voro::container &con, std::vector<int> cells, const double &alfa, const double &beta, const double &B);
bool feasibility(voro::container_poly &con, std::vector<int> cells, const double &alfa, const double &beta, const double &B);
// ^ -- for a given list of cells verifies its feasibility w.r.t hardcore parameters (feasibility is the property of each cell)

bool feasibility(voro::container &con, int ijk1, int q1, const double &alfa, const double &beta, const double &B);   
// ^ -- feasibility function for special case ADD (verifies feasibility for a neighborhood of added particle after adding)
bool feasibility(voro::container &con, std::vector<int> &rsr, std::vector<int> &rsio, std::vector<double> &rsap, const double &alfa, const double &beta, const double &B); 
// ^ -- feasibility function for special case DELETE (verifies feasibility for a neighborhood of deleted particle after deleting)

void bdma_step(long &npart, std::vector<int> &fid, voro::container &con, double sigma, double theta, double zet, 
	const double &alfa, const double &beta, const double &B, long &rn_add, long &rn_del, long &rn_mov);
// ^ -- the main body of the BDM algorithm; with probability 1/3 tries to add/delete or move a generator, examines the feasibility of the operation
// and computes the energy of such operation (this is done by calling appropriate try_function - section try functions below)
// if the operation is accepted, the function modifies the container permanently
void bdma_step_2(long &npart, std::vector<int> &fid, voro::container &con, double sigma, double theta, double zet, 
	const double &alfa, const double &beta, const double &B, long &rn_add, long &rn_del, long &rn_mov, bool &fc);
void bdma_step_3(long &npart, std::vector<int> &fid, voro::container &con, double sigma, double theta1, double theta2, double zet,
	const double &alfa, const double &beta, const double &B, long &rn_add, long &rn_del, long &rn_mov, bool &fc);



double V1(voro::container &con, int pr, std::vector<int> &rsr);
// ^ -- computes mean number of faces for the specified region of cells; designed for energy modifying the number of faces of the average cell
double V1_nf(voro::voronoicell_neighbor &rc, int &n0);
// ^ -- energy preferring such cells which number of neighbours is similar to the observed mean number of neighbours
double V2(voro::voronoicell_neighbor &rc1, voro::voronoicell_neighbor &rc2, double &rx, double &ry, double &rz, double &rxx, double &ryy, double &rzz);
// ^ -- for two cells C1, C2 sharing a face returns sqrt(max(vol C1, vol C2) / min(vol C1, vol C2) - 1) iff the barycenter of the union of C1 and C2 lies inside 
// the window (this is tested by calling the function barycentrum; (vol = volume)
double V3f(int &id1, voro::voronoicell_neighbor &rc1, int &id2, voro::voronoicell_neighbor &rc2, int &id3, voro::voronoicell_neighbor &rc3,
	double &rx, double &ry, double &rz, double &rxx, double &ryy, double &rzz, double &rxxx, double &ryyy, double &rzzz);
// ^ -- for three cells C1, C2, C3 sharing an edge returns sqrt(max(far f1, far f2, far f3) / min (far f1, far f2, far f3) - 1) iff the barycenter of the
// union C1, C2 and C3 lies inside the window; (far = area of the face containing the shared edge (there are just three such faces))
double V3a(int &id1, voro::voronoicell_neighbor &rc1, int &id2, voro::voronoicell_neighbor &rc2, int &id3, voro::voronoicell_neighbor &rc3,
	double &rx, double &ry, double &rz, double &rxx, double &ryy, double &rzz, double &rxxx, double &ryyy, double &rzzz);
// ^ -- for three cells C1, C2, C3 sharing an edge returns sqrt(max(dan C1, dan C2, dan C3) / min (dan C1, dan C2, dan C3) - 1) iff the barycenter of the
// union C1, C2 and C3 lies inside the window; (dan = the size of the dihedral angle belonging to the shared edge)
double V4(int &id1, voro::voronoicell_neighbor &rc1, int &id2, voro::voronoicell_neighbor &rc2, int &id3, int &id4);

void bdma_step(long &npart, std::vector<int> &fid, voro::container_poly &conp, double sigma, double theta, double zet, const double &alfa, const double &beta, 
	const double &B, const double &iota, long &rn_add, long &rn_del, long &rn_mov, double &arad);
// ^ -- BDM algorithm adapted to the Laguerre case; besides adding/deleting or moving the generator (i.e. its coordinates) enables the change of the 
// radius, each of these operations is done with probability 1/4



// try functions //
// - returns the exp of energy necessary to do appropriate operation; firstly verifies the feasibility of the operation and if feasible computes the energy;
// returns 0 iff not feasible; both feasibility and the change of energy is computed locally; following fcs differs by used energy function and for each type
// there are two versions (one for Voronoi tess. and another for Laguerre)
// try functions I - use V2 as the energy function
double try_add(int id, double &nx, double &ny, double &nz, voro::container &con, const double &theta, const double &alfa, const double &beta, const double &B);
double try_delete(int &ijk_del, int &q_del, voro::container &con, const double &theta, const double &alfa, const double &beta, const double &B);
double try_move(int &ijk_del, int &q_del, double &nx, double &ny, double &nz, voro::container &con, const double &theta, const double &alfa, const double &beta, const double &B);

double try_add(int id, double &nx, double &ny, double &nz, double &nrad, voro::container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, const double &iota);
double try_delete(int &ijk_del, int &q_del, voro::container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, const double &iota);
double try_MOVE(int &ijk_del, int &q_del, double &nx, double &ny, double &nz, double &nrad, voro::container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, const double &iota);
double try_move(int &ijk_del, int &q_del, double &nx, double &ny, double &nz, voro::container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, const double &iota);
double try_change_rad(int &ijk_del, int &q_del, double &nrad, voro::container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, const double &iota);

// try functions II - use V3f/V3a as the energy function
double try_add_2(int id, double &nx, double &ny, double &nz, voro::container &con, const double &theta, const double &alfa, const double &beta, const double &B, bool &fc);
double try_delete_2(int &ijk_del, int &q_del, voro::container &con, const double &theta, const double &alfa, const double &beta, const double &B, bool &fc);
double try_move_2(int &ijk_del, int &q_del, double &nx, double &ny, double &nz, voro::container &con, const double &theta, const double &alfa, const double &beta, const double &B, bool &fc);

double try_add_2(int id, double &nx, double &ny, double &nz, double &nrad, voro::container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, bool &fc);
double try_delete_2(int &ijk_del, int &q_del, voro::container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, bool &fc);
double try_move_2(int &ijk_del, int &q_del, double &nx, double &ny, double &nz, voro::container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, bool &fc);
double try_change_rad_2(int &ijk_del, int &q_del, double &nrad, voro::container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, bool &fc);

// try functions III - goal should be to obtain energy preferring configurations with given mean number of cell faces
double try_add_e(int nf, int id, double &nx, double &ny, double &nz, voro::container &con, const double &theta, const double &alfa, const double &beta, const double &B);
double try_delete_e(int nf, int &ijk_del, int &q_del, voro::container &con, const double &theta, const double &alfa, const double &beta, const double &B);
double try_move_e(int nf, int &ijk_del, int &q_del, double &nx, double &ny, double &nz, voro::container &con, const double &theta, const double &alfa, const double &beta, const double &B);



// statistics //
double corr_coef(voro::container &rcon);
double corr_coef(voro::container_poly &rcon);
// ^ -- returns the correlation coefficient of the tessellation
void con_stats(voro::container_poly &con);
void cell_stats(voro::container &rcon);
void cell_stats(voro::container_poly &rcon);
// ^ -- produces the output into cell_stats.txt of cell statistics such as volume, surface area, total edge distance, number of faces, number of edges,
// number of vertices, second power of the maximal distance between the face and generator
void neigh_list(voro::container_poly &rcon);
// ^ -- outputs the neighbors list for each particle
void face_stats(voro::container &rcon);
void face_stats(voro::container_poly &rcon);
// ^ -- produces the output into face_stats.txt of face statistics such as face area, number of edges/vertices, perimeter, volumes of the two cells sharing this face
void edge_stats(voro::container &rcon);
void edge_stats(voro::container_poly &rcon);
// ^ -- produces the output into edge_stats.txt of edge statistics such as length, areas of three faces sharing this edge, three dihedral angles (in radians)
void vertex_stats(voro::container &rcon);
void vertex_stats(voro::container_poly &rcon);
void stats_output(voro::container &rcon);

void face_edge_lengths(std::vector<double> &rv, voro::voronoicell_neighbor &rc);
void edge_lengths(std::vector<double> &rv, voro::voronoicell_neighbor &rc);
void dihedral_angles(std::vector<double> &rv, voro::voronoicell_neighbor &rc);



// estimation //
void hardcore_estim(voro::container &rcon, double &ralfa_e, double &rbeta_e, double &rB_e);
void hardcore_estim(voro::container_poly &rcon, double &ralfa_e, double &rbeta_e, double &rB_e);
// ^ -- returns alfa_e, beta_e, B_e - the estimates of hardcore parameters
void overlap_e(voro::container_poly &rcon, double &riota_e);

int no_removable(voro::container &rcon, const double &ralfa, const double &rbeta, const double &rB, std::vector<double> &vb, double lx, double ux, double ly, double uy, double lz, double uz);
int no_removable(voro::container_poly &rcon, const double &ralfa, const double &rbeta, const double &rB, const double &riota, std::vector<double> &vb, double lx, double ux, double ly, double uy, double lz, double uz);
// ^ -- returns the number of removable points (i.e. generators) within the specified subwindow (bounded by lb and ub in each coordinate)
// moreover computes the energy necessary to delete removable points (these energies are stored in vector vb)
int no_removable_2(voro::container &rcon, const double &ralfa, const double &rbeta, const double &rB);
int no_removable_2(voro::container_poly &rcon, const double &ralfa, const double &rbeta, const double &rB);
// ^ -- another way how to compute the number of removable points in the whole container (i.e. whole window)
bool removable(voro::container &lcon, int rid, const double &ralfa, const double &rbeta, const double &rB); // pouzita v residuich
bool removable(voro::container_poly &lcon, int rid, const double &ralfa, const double &rbeta, const double &rB);
bool addable(voro::container &rcon, double &rx, double &ry, double &rz, const double &ralfa, const double &rbeta, const double &rB);
bool addable(voro::container_poly &rcon, double &rx, double &ry, double &rz, double &rr, const double &ralfa, const double &rbeta, const double &rB, const double &riota);
// ^ -- returns true if the given point can be added into container (i.e. if the tessellation with a new generator is still feasible)
// calls feasibility function for a special case ADD
int no_points(voro::container &rcon, double lb, double ub);
int no_points(voro::container_poly &rcon, double lb, double ub);
// ^ -- returns number of points in the specified subwindow
int total_points(voro::container &rcon, const double &ralfa, const double &rbeta, const double &rB); // kontrolni, nepouzita

void estim(double &th_estim, double &z_estim, voro::container &rcon, int rN, double &ralfa_e, double &rbeta_e, double &rB_e, double lx, double ux, double ly, double uy, double lz, double uz);
// ^ -- the main function of the estimation section; returns the estimates th_estim and z_estim of parameters theta and zet; works on previously specified
// subwindow; calls fc no_removable (determines number of removable points and its energies needed to deletion), then generates random points and computes 
// its local energies (calls addable to determine if the point can be added), these points are used for approximating integrals by MC method;
// precomputed energies (stored in vectors) are then used in calling of the Newton-Raphson procedure, where specify the function whose root are searched;
// this root is estimate of theta and is used afterwards in determining estimate of zet; these estimates can be founded in many ways, here implemented are
// so called "quick" estimate (quick does not reffer to the time consumed during estimation) and max pseudolikelihood estimate (this can be obtained in 
// two ways - first by Newton-Raphson searching the root of derivative of PLL function; or by LBFGS optimalization method searching directly min of PLL

// note: this part is still under development and therefore changing rapidly
void estim(double &th_estim, double &z_estim, voro::container_poly &rcon, int N, double &ralfa_e, double &rbeta_e, double &rB_e, double &riota_e, double lx, double ux, double ly, double uy, double lz, double uz);


double zet_estim(int rN, double &rth, std::vector<bool> &ra, std::vector<double> &rc); // nepouzita
double PLL_const(voro::container &rcon, double &ralfa_e, double &rbeta_e, double &rB_e); // nepouzita
long PLL_coeff(voro::container &rcon, double &ralfa_e, double &rbeta_e, double &rB_e, int rN, std::vector<bool> &ra, std::vector<double> &rc); // nepouzita

void raw_resid(std::vector<double> &raw_res, voro::container &rcon, int &gsi, int rN, double &ralfa_e, double &rbeta_e, double &rB_e, double &th_estim, double &z_estim);
void raw_resid(std::vector<double> &raw_res, voro::container_poly &rcon, int &gsi, int rN, double &ralfa_e, double &rbeta_e, double &rB_e, double riota_e, double &th_estim, double &z_estim);
void inverse_resid(std::vector<double> &inv_res, voro::container &rcon, int &gsi, int rN, double &ralfa_e, double &rbeta_e, double &rB_e, double &th_estim, double &z_estim);
void pearson_resid(std::vector<double> &p_res, voro::container &rcon, int &gsi, int rN, double &ralfa_e, double &rbeta_e, double &rB_e, double &th_estim, double &z_estim);
void smooth_res_field(voro::container &rcon, double &ralfa_e, double &rbeta_e, double &rB_e, double &th_estim, double &z_estim);
// ^ -- fcs computing residual fields of tessellation

void MLE(double th_estim, double alfa_e, double beta_e, double B_e, double zet_e);
// ^ ve vyvoji

double con_energy(voro::container &rcon, double lb, double ub);
double con_energy(voro::container_poly &rcon, double lb, double ub);





//
//

//
//

class histogram {
public: // beginning value = the smallest value to be stored
	double sp;
	// step size
	double step;
	// vector of occurencies (for a given range);
	std::vector<int> oc;
	// number of occurencies for modus 
	double ocm;
	// sum of all occurencies
	double so;
	// number of occurencies not included in oc
	double noc;  // allows to compare two histograms more easily


	histogram() { sp = 0; step = 1; oc.push_back(1); ocm = 1; so = 1; noc = 0; };
	histogram(double, double, std::vector<int>);

	void read_hist_vol();
	void read_hist_nof();
	void create_hist_int(voro::container_poly &con);
	void create_hist_double(voro::container_poly &con);
	int hist_value(double rc);
	void hist_act(double val, bool op);
};

// void hist_act(histogram &hist, double val, bool op);
double hist_dis(histogram &hist1, histogram &hist2);
double hist_disp(histogram &hist1, histogram &hist2);


class con_info {
public: 
	// total particles  -- prtz se mohou vyskytovat prazdne bunky, musi se ukladat informace pred i po
	int tp_bef, tp_aft;
	// sample mean (possibly more sample means --> vector) (the sum)
	//double mean_bef; // 
	std::vector<double> mean_bef;
	//double mean_aft; // 
	std::vector<double> mean_aft;
	// sample variance (the sum)
	std::vector<double> var_bef;
	//double mean_aft; // 
	std::vector<double> var_aft;
	// ...
	// histogram (possibly histograms)
	histogram hist_bef;
	histogram hist_aft;
	histogram hist2_bef;
	histogram hist2_aft;
	// empty cells
	std::vector<int> empty;

	con_info() { tp_bef = 0; tp_aft = 0; mean_bef.push_back(0); mean_aft.push_back(0); var_bef.push_back(0); var_aft.push_back(0); empty.clear(); };
	con_info(int, std::vector<double>, std::vector<double>, histogram);
	void clear() { mean_bef.clear(); mean_aft.clear(); var_bef.clear(); var_aft.clear(); };

	// !!!dodrzuj poradi: 
	void get_tp(voro::container_poly &con) { tp_bef = nonempty_cells(con); tp_aft = tp_bef; };
	//void get_mean(voro::container_poly &con);
	void get_meansum(voro::container_poly &con);
	//void get_hist --> use functions create hist
	void get_varsum(voro::container_poly &con);
	void varsum(voro::container_poly &con);

};

void delete_empty(voro::container_poly &con);



// LAGUERRE - snaha o efektivni algoritmus

void LAG_bdma_cout(long &npart, std::vector<int> &fid, voro::container_poly &conp, voro::container_poly &conp_copy, double sigma, std::vector<double> theta, double zet, std::vector<double> hard_par, long &rn_add, long &rn_del, long &rn_mov, con_info &info);
void LAG_bdma(long &npart, std::vector<int> &fid, voro::container_poly &conp, voro::container_poly &conp_copy, double sigma, std::vector<double> theta, double zet, std::vector<double> hard_par, long &rn_add, long &rn_del, long &rn_mov, con_info &info, histogram &hist, histogram &hist2);

double V_function(double val);
double V_function(double val, histogram &hist);
double V_function_2(double val, histogram &hist);
double V_function(histogram &val, histogram &hist);
double V_function_2(histogram &val, histogram &hist);
double V_function(double val1, double val2);
double V_function(double val1, double val2, histogram &hist);

double V2_function(voro::voronoicell_neighbor &rc1, voro::voronoicell_neighbor &rc2);

double LAG_recompute(voro::container_poly &con, voro::container_poly &newcon, int type, int id, std::vector<double> &h_par, std::vector<double> &theta, con_info &info);

void LAG_container(voro::container_poly &con, voro::container_poly &newcon, int type, int id);
bool LAG_cells(voro::container_poly &con, voro::container_poly &newcon, int type, int id, con_info &info, std::vector<int> &cells, std::vector<int> &cells_pos);
void LAG_sec(voro::container_poly &con, voro::container_poly &newcon, int id, std::vector<int> cells, std::vector<int> cells_pos, std::vector<int> &sec, std::vector<int> &sec_pos);
bool LAG_feasibility(voro::container_poly &con, voro::container_poly &newcon, std::vector<double> &h_par, std::vector<int> cells_pos);
void LAG_V1(voro::container_poly &con, voro::container_poly &newcon, int type, int id, con_info &info, histogram &hist, histogram &hist2, std::vector<double> &parts, std::vector<int> cells, std::vector<int> cells_pos);
double LAG_V2(voro::container_poly &con, voro::container_poly &newcon, int type, int id, con_info &info, std::vector<int> cells, std::vector<int> cells_pos, std::vector<int> sec, std::vector<int> sec_pos);

void LAG_estim(double &th_estim, double &z_estim, voro::container_poly &con, voro::container_poly &newcon, int N, std::vector<double> &hpar_estim, con_info &info, histogram &hist, histogram &hist2, double lx, double ux, double ly, double uy, double lz, double uz);
int LAG_removable(voro::container_poly &con, voro::container_poly &con_copy, std::vector<double> &vb, std::vector<double> &hpar_estim, con_info &info, histogram &hist, histogram &hist2, double lx, double ux, double ly, double uy, double lz, double uz);

// include LAG_recompute, dens_ener, ...

