#include "Header.h"

using namespace voro;

/*// example of class-function
class Rosenbrock
{
private:
	int n;
public:
	Rosenbrock(int n_) : n(n_) {}
	double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
	{
		double fx = 0.0;
		for (int i = 0; i < n; i += 2)
		{
			double t1 = 1.0 - x[i];
			double t2 = 10 * (x[i + 1] - x[i] * x[i]);
			grad[i + 1] = 20 * t2;
			grad[i] = -2.0 * (x[i] * grad[i + 1] + t1);
			fx += t1 * t1 + t2 * t2;
		}
		return fx;
	}
};

// class PLL (pseudologlikelihood function) allow us to create function PLL suitable for LBFGS optimalization algorithm
class PLL
{
private:
	int n;
	std::vector<int> a;
	std::vector<double> b, c;
	int k;
public:
	// constructor takes values of function coefficients
	PLL(int n_, std::vector<int> a_, std::vector<double> b_, std::vector<double> c_, int k_)
	{
		n = n_; a = a_; b = b_; c = c_; k = k_;
	}
	// definition of the PLL function and computation of its gradient
	double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
	{
		double fx = 0.0;
		int i, j;
		double t1 = 0;

		for (j = 0; j < k; j++) { t1 = t1 + b[j]; }

		for (i = 0; i < n; i += 2)
		{
			double t2 = 0;
			double t3 = 0;
			for (j = 0; j < c.size(); j++) { 
				t2 = t2 + exp(-x[i + 1] * c[j])*a[j]; 
				t3 = t3 + c[j] * exp(-x[i + 1] * c[j])*a[j];
			}

			//double t1 = 1.0 - x[i];
			//double t2 = 10 * (x[i + 1] - x[i] * x[i]);
			grad[i + 1] = -x[i] * t3 + t1;
			grad[i] = t2 - k / x[i];
			//grad[i + 1] = 20 * t2;
			//grad[i] = -2.0 * (x[i] * grad[i + 1] + t1);
			fx += x[i] * t2 + x[i + 1] * t1 - log(x[i])*k;
			//fx += t1 * t1 + t2 * t2;
		}
		return fx;
	}
};
*/


void hardcore_estim(voro::container &rcon, double &ralfa_e, double &rbeta_e, double &rB_e)
{
	// [in]		con				the container with stored particles.
	// [out]	alfa_e			the estimate of hardcore parameter alfa
	// [out]	beta_e			the estimate of hardcore parameter beta
	// [out]	B_e				the estimate of hardcore parameter B

	unsigned int k, fng;
	int id, i, j;
	double x, y, z, xn, yn, zn;
	double dist, vol;
	voronoicell_neighbor c;  // bunka
	std::vector<int> neigh;  // vektor pro prirazeni ID sousedu
	std::vector<int> vert;	 // vektor obsahujici vrcholy kazde steny bunky v nejakem danem poradi

	ralfa_e = 2;
	rbeta_e = 0;
	rB_e = 0;

	// pro prochazeni pouziji tridu loop
	//c_loop_all clo(rcon);          // loop pres vsechny castice
								   // c_loop_order clo(con, po);    // loop jen pres vybrane

	//if (clo.start()) do if (rcon.compute_cell(c, clo)) {
		// start() - sets the class to consider the first particle 
		// compute_cell(c, clo) - return true if the computation of the Voronoi cell for a particle currently being referenced by a loop class was succesful 

		//id = clo.pid();     // Get the ID of the current particle under consideration
		//clo.pos(x, y, z);   // Get the position of the current particle under consideration

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			id = rcon.id[j][i];
			x = rcon.p[j][3 * i];
			y = rcon.p[j][3 * i + 1];
			z = rcon.p[j][3 * i + 2];

			rcon.compute_cell(c, j, i);

			c.neighbors(neigh);  // Computes a vector list of neighbors of the current particle under consideration
			c.face_vertices(vert);					// {PER} computes list of vertices for each face of the cell
			vol = c.volume();

			fng = 0;

			// loop over the neighbors
			for (k = 0; k < neigh.size(); k++) {

				// Skip if the neighbor information is smaller than this particle's ID, to avoid double counting. This also
				// removes faces that touch the walls, since the neighbor information is set to negative numbers for these cases.
				if (neigh[k] > id) {

					// vzdalenost dvou bodu: sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 ) --> fce point_dist
					face_dist(fng, vert, x, y, z, xn, yn, zn, c);		// urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
					xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;	// urcim skutecne souradnice tohoto souseda za stenou

					dist = point_dist(x, y, z, xn, yn, zn);
					dist = dist / 2;
					if (dist < ralfa_e) { ralfa_e = dist; }   // actualization of estimates
					if (dist > rbeta_e) { rbeta_e = dist; }
					dist = pow(dist, 3) / vol;
					if (dist > rB_e) { rB_e = dist; }
				}

				fng = fng + vert[fng] + 1;			// set actual position in vector of face vertices
			}
		}
	}
	//} while (clo.inc());    // inc() - finds the next particle to test.
}

void hardcore_estim(voro::container_poly &rcon, double &ralfa_e, double &rbeta_e, double &rB_e)
{
	// [in]		con				the container with stored particles.
	// [out]	alfa_e			the estimate of hardcore parameter alfa
	// [out]	beta_e			the estimate of hardcore parameter beta
	// [out]	B_e				the estimate of hardcore parameter B

	int id, id2, i, j, k, l;
	double x, y, z, xx, yy, zz, tx, ty, tz, r, rr;
	double h_min, h_max, vol, dist;
	bool grain;
	voronoicell_neighbor c;  // bunka
	
	ralfa_e = 200000;
	rbeta_e = 0;
	rB_e = 0;

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			id = rcon.id[j][i];
			x = rcon.p[j][4 * i];
			y = rcon.p[j][4 * i + 1];
			z = rcon.p[j][4 * i + 2];

			grain = rcon.compute_cell(c, j, i);
			if (grain == true) {		// nonempty cell
				c.centroid(tx, ty, tz);					// centroid in relative coordinates wrt the generator
				tx = x + tx; ty = y + ty; tz = z + tz;	// centroid in absolute coordinates

				vol = c.volume();

				h_fcs(c, x, y, z, tx, ty, tz, h_max, h_min);

				if (h_min < ralfa_e) { ralfa_e = h_min; }   // actualization of estimates
				if (h_max > rbeta_e) { rbeta_e = h_max; }
				h_max = pow(h_max, 3) / vol;
				if (h_max > rB_e) { rB_e = h_max; }

			}	// END..if(nonempty cell)
		}
	}
}


void overlap_e(voro::container_poly &rcon, double &riota_e)
{
	// [in]		con				the container with stored particles.
	// [out]	riota_e			the estimate of hardcore parameter iota
	
	int i, j, k, l, id1, id2;
	double dist, x, y, z, xx, yy, zz, tx, ty, tz, r, rr;

	riota_e = 200000;

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			id1 = rcon.id[j][i];

			x = rcon.p[j][4 * i];
			y = rcon.p[j][4 * i + 1];
			z = rcon.p[j][4 * i + 2];

			r = rcon.p[j][4 * i + 3];

			for (l = 0; l < rcon.nxyz; l++) { // loop over boxes
				for (k = 0; k < rcon.co[l]; k++) { // loop over particles in considered box

					id2 = rcon.id[l][k];

					if (id1 < id2) {					// prevent doublecounting and using one generator two-times
						xx = rcon.p[l][4 * k];
						yy = rcon.p[l][4 * k + 1];
						zz = rcon.p[l][4 * k + 2];

						dist = point_dist(x, y, z, xx, yy, zz);

						if (dist == 0 ){std::cout << x << " " << y << " " << z << " ; " << xx << " " << yy << " " << zz << " ; " << id1 << " " << id2 << "\n";}

						rr = rcon.p[l][4 * k + 3];

						if (riota_e > (dist - r)) { riota_e = dist - r; }
						if (riota_e > (dist - rr)) { riota_e = dist - rr; }
					}
				}
			} // end second loop		

		}
	} // end first loop
}

// fc no_removable returns number of removable points in the container (i.e. number of points from container which
//   when removed the new configuration is still feasible
int no_removable(voro::container &rcon, const double &ralfa, const double &rbeta, const double &rB, std::vector<double> &vb, double lx, double ux, double ly, double uy, double lz, double uz)
{
	// [in]		con				the container with stored particles
	// [in]     alfa, beta, B	feasibility constraints
	// [out]	rb				vector of coefficients
	// [in]		lb, ub			the lower and upper bound of considered subwindow

	int id, no, i, j, ijk, q, ci;
	double eloc = 0;

	//	c_loop_all clo(rcon);          // loop pres vsechny castice
	// loop se neumi vyporadat s tim, ze fc removable meni strukturu containeru
	no = 0; 

	/* pro kazdy bod chceme nyni overit zda konfigurace po jeho odebrani je pripustna, mozno dvema zpusoby:
	1/ fce feasibility pro specialni pripad smazani bodu - po smazani bodu zavolam tuto specialni fci, abych tak mohl
	ucinit, musim ale nejdriv spocitat vektory sr, sio, sap
	bool feasibility(voro::container &con, std::vector<int> &rsr, std::vector<int> &rsio, std::vector<double> &rsap, const double &alfa, const double &beta)
	2/ fce feasibility pro pocatecni konfiguraci - po smazani bodu overim celou konfiguraci; ne moc efektivni, nebot
	vetsina bunek odstranenim jedne z nich nebyla ovlivnena
	bool feasibility(container &con, const double &alfa, const double &beta) // omezeni - mely by se predat fci odkazem jako konstanty
	*/
	// fc removable vyuzije prvni pristup

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		ci = 0; // poradove cislo castice
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box - prtz se castice smaze a pak zase prida,
										   // je potreba co[j]-krat vzit prvni castici v bodu, ta se vzdy presune na konec
			id = rcon.id[j][ci];
			// std::cout << rcon.co[j] << " "; ///////////////////////////////////////////////////////////////////////////
			//std::cout << id << " ";  /////////////////////////////////////////////////////////////////////////////////
			//if (removable(rcon, id, ralfa, rbeta, rB)) { 
				find_pos(ijk, q, id, &rcon);	// find position of the particle

				if (rcon.p[ijk][3*q] > lx && rcon.p[ijk][3*q] < ux && rcon.p[ijk][3*q + 1] > ly && rcon.p[ijk][3*q + 1] < uy && rcon.p[ijk][3*q + 2] > lz && rcon.p[ijk][3*q + 2] < uz) {
					// if the coordinates are in subwindow

					eloc = try_delete(ijk, q, rcon, 1, ralfa, rbeta, rB);
					if (eloc > 0) {
						vb.push_back(log(eloc)); // save this information
						no++;
						//std::cout << id << " "; /////////////////////////////////////////////////////////////////////////////////					
					}  
				} else { ci++; }
			//ci++;
			//}
		}		
	}
	//std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////
	//std::cout << ci << "\n";
	return no;
}

int no_removable(voro::container_poly &rcon, const double &ralfa, const double &rbeta, const double &rB, const double &riota, std::vector<double> &vb, double lx, double ux, double ly, double uy, double lz, double uz)
{
	// [in]		con				the container with stored particles
	// [in]     alfa, beta, B	feasibility constraints
	// [out]	vb				vector of coefficients
	// [in]		lb, ub			the lower and upper bound of considered subwindow (in every coordinate)

	int id, no, i, j, ijk, q, ci;
	double eloc = 0;

	no = 0;

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		ci = 0; // poradove cislo castice; indikator poctu neodstranitelnych bodu = pouze pro kontrolu, zda-li jde o doplnek do celkoveho poctu
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box - prtz se castice smaze a pak zase prida,
										   // je potreba co[j]-krat vzit prvni castici v bodu, ta se vzdy presune na konec
			id = rcon.id[j][ci];
			// std::cout << rcon.co[j] << " "; ///////////////////////////////////////////////////////////////////////////
			//std::cout << id << " ";  /////////////////////////////////////////////////////////////////////////////////
			//if (removable(rcon, id, ralfa, rbeta, rB)) { 
			find_pos(ijk, q, id, &rcon);	// find position of the particle

			if (rcon.p[ijk][4 * q] > lx && rcon.p[ijk][4 * q] < ux && rcon.p[ijk][4 * q + 1] > ly && rcon.p[ijk][4 * q + 1] < uy && rcon.p[ijk][4 * q + 2] > lz && rcon.p[ijk][4 * q + 2] < uz) {
				// if the coordinates are in subwindow

				eloc = try_delete(ijk, q, rcon, 1, ralfa, rbeta, rB, riota);
				if (eloc > 0) {
					vb.push_back(log(eloc)); // save this information
					no++;
					//std::cout << id << " "; /////////////////////////////////////////////////////////////////////////////////
				}
			}
			else { ci++; }
			//ci++;
			//}0
		}
	}
	//std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////
	//std::cout << ci << "\n";
	return no;
}

// mene efektivni fce overujici pripustnost cele mozaiky
int no_removable_2(voro::container &rcon, const double &ralfa, const double &rbeta, const double &rB)
{
	// [in]		con				the container with stored particles
	// [in]     alfa, beta, B	feasibility constraints
	// [in]		lb, rb			the lower and upper bound of considered subwindow

	int id, no, i, j;
	int nul;
	double x, y, z;
	no = 0;

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box - prtz se castice smaze a pak zase prida,
										   // je potreba co[j]-krat vzit prvni castici v bodu, ta se vzdy presune na konec
			id = rcon.id[j][0];
			x = rcon.p[j][0];			// computes coordinates of the particle
			y = rcon.p[j][1];
			z = rcon.p[j][2];
			// std::cout << rcon.co[j] << " "; ///////////////////////////////////////////////////////////////////////////
			// std::cout << id << "\n";  /////////////////////////////////////////////////////////////////////////////////
			nul = 0;
			erase(j, nul, &rcon);
			if (feasibility(rcon, ralfa, rbeta, rB)) { no++; } // pozor vypisuje
			rcon.put(id, x, y, z);	// navraceni castice
		}
	}

	return no;
}

int no_removable_2(voro::container_poly &rcon, const double &ralfa, const double &rbeta, const double &rB)
{
	// [in]		con				the container with stored particles
	// [in]     alfa, beta, B	feasibility constraints
	// [in]		lb, rb			the lower and upper bound of considered subwindow

	int id, no, i, j;
	int nul;
	double x, y, z, r;
	no = 0;

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box - prtz se castice smaze a pak zase prida,
										   // je potreba co[j]-krat vzit prvni castici v bodu, ta se vzdy presune na konec
			id = rcon.id[j][0];
			x = rcon.p[j][0];			// computes coordinates of the particle
			y = rcon.p[j][1];
			z = rcon.p[j][2];
			r = rcon.p[j][3];
			// std::cout << rcon.co[j] << " "; ///////////////////////////////////////////////////////////////////////////
			// std::cout << id << "\n";  /////////////////////////////////////////////////////////////////////////////////
			nul = 0;
			erase(j, nul, &rcon);
			if (feasibility(rcon, ralfa, rbeta, rB)) { no++; } // pozor vypisuje
			rcon.put(id, x, y, z, r);	// navraceni castice
		}
	}

	return no;
}

// fc removable returns true if for the given point is possible to be removed without breaking the feasibility condition
bool removable(voro::container &lcon, int rid, const double &ralfa, const double &rbeta, const double &rB)
{
	// [in]		con				the container with stored particles
	// [in]		id				ID of the point under examination
	// [in]     alfa, beta, B	feasibility constraints

	unsigned int fng, k, i;
	int ijk_del, q_del, ijk, q;
	double x, y, z, xn, yn, zn;
	voronoicell_neighbor c;
	std::vector<int> neigh;		// vektor pro prirazeni ID sousedu
	std::vector<int> vert;
	std::vector<int> sr;		// vectors storing ijk,q information about secondary particles
	std::vector<int> sio;		// information in/out for secondary particles
	std::vector<double> sap;	// real positions of secondary particles which are "out"

	find_pos(ijk_del, q_del, rid, &lcon);	// find position of the particle
	lcon.compute_cell(c, ijk_del, q_del);

	x = lcon.p[ijk_del][3 * q_del];			// computes coordinates of the particle
	y = lcon.p[ijk_del][3 * q_del + 1];
	z = lcon.p[ijk_del][3 * q_del + 2];

	c.neighbors(neigh);						// compute its neighbors			
	c.face_vertices(vert);					// computes list of vertices for each face of the cell
	fng = 0;								//
	k = 0;

	// now we would like to save structural vectors sr, sio, sap
	for (i = 0; i < neigh.size(); i++) {	// loop over the neighbors
		find_pos(ijk, q, neigh[i], &lcon);  // find position of the i-th neighbor

		face_dist(fng, vert, x, y, z, xn, yn, zn, c);		// urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
		xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;	// urcim skutecne souradnice tohoto souseda za stenou

		sr.push_back(ijk); sr.push_back(q); // add this neighbor to the second particle vector

		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		 // je-li pozice souseda v okne (tj. nedoslo k preklopeni)
			sio.push_back(0);
		}
		else {
			k++;
			sio.push_back(k);					// neighbor is k-th "out"
			sap.push_back(xn); sap.push_back(yn); sap.push_back(zn);  // storing real position of this neighbor
		}

		fng = fng + vert[fng] + 1;			// set actual position in vector of face vertices
	}

	erase(ijk_del, q_del, &lcon);

	if (feasibility(lcon, sr, sio, sap, ralfa, rbeta, rB)) { lcon.put(rid, x, y, z); return 1; }
	else { lcon.put(rid, x, y, z); return 0; }
}

bool removable(voro::container_poly &lcon, int rid, const double &ralfa, const double &rbeta, const double &rB)
{
	// [in]		con				the container with stored particles
	// [in]		id				ID of the point under examination
	// [in]     alfa, beta, B	feasibility constraints

	unsigned int fng, k, i;
	int ijk_del, q_del, ijk, q;
	double x, y, z, r, xn, yn, zn;
	voronoicell_neighbor c;
	std::vector<int> neigh;		// vektor pro prirazeni ID sousedu
	std::vector<int> vert;
	std::vector<int> sr;		// vectors storing ijk,q information about secondary particles
	std::vector<int> sio;		// information in/out for secondary particles
	std::vector<double> sap;	// real positions of secondary particles which are "out"

	find_pos(ijk_del, q_del, rid, &lcon);	// find position of the particle
	lcon.compute_cell(c, ijk_del, q_del);

	x = lcon.p[ijk_del][4 * q_del];			// computes coordinates of the particle
	y = lcon.p[ijk_del][4 * q_del + 1];
	z = lcon.p[ijk_del][4 * q_del + 2];
	r = lcon.p[ijk_del][4 * q_del + 3];

	erase(ijk_del, q_del, &lcon);

	if (feasibility(lcon, ralfa, rbeta, rB)) { lcon.put(rid, x, y, z, r); return 1; }
	else { lcon.put(rid, x, y, z, r); return 0; }
}

// fc addable returns true if for the given point is possible to be added without breaking the feasibility condition
// - mene efektivni fce overujici pripustnost cele mozaiky (tak jako fce no_removable_2)
bool addable(voro::container &rcon, double &rx, double &ry, double &rz, const double &ralfa, const double &rbeta, const double &rB)
{
	// [in]		con				the container with stored particles
	// [in]		x,y,z			point under examination
	// [in]     alfa, beta, B	feasibility constraints - hardcore parameters

	int i;
	int ijk_add, q_add, ijk, q;
	std::vector<int> cells;
	std::vector<int> neigh;
	voronoicell_neighbor c;

	// chci pridat bod, k tomu potrebuji nejake ID, jak ale zjistim, ktere ID je dostupne? co kdybych pouzil jako ID nejake
	//   cislo, ktere se bezne jako ID nepouziva (treba 0 nebo neco zaporneho)

	rcon.put(0, rx, ry, rz);			// add particle under examination
	find_pos(ijk_add, q_add, 0, &rcon); // find position of added particle

	rcon.compute_cell(c, ijk_add, q_add);
	c.neighbors(neigh);

//  1) Feasibility for the neighborhood pr+sr:
//	cells.clear();

//	for (i = 0; i < neigh.size(); i++) {
//		find_pos(ijk, q, neigh[i], &rcon);
//		cells.push_back(ijk); cells.push_back(q); // add this neighbor to the vector
//	}
//	cells.push_back(ijk_add); cells.push_back(q_add);

//	if (feasibility(rcon, cells, ralfa, rbeta, rB)) { erase(ijk_add, q_add, &rcon); return 1; }

//  2) Feasibility for the whole container:
	if (feasibility(rcon, ralfa, rbeta, rB)) { erase(ijk_add, q_add, &rcon); return 1; }

	// verify feasibility:
	//if (feasibility(rcon, ralfa, rbeta, rB)) { erase(ijk_add, q_add, &rcon);  return 1; }  // nemel by se tu pouzit spec pripad feasibility pro add
	//if (feasibility(rcon, ijk_add, q_add, ralfa, rbeta, rB)) { erase(ijk_add, q_add, &rcon);  return 1; }  // pouzita feasibility pro ADD

	erase(ijk_add, q_add, &rcon);
	return 0;
}


bool addable(voro::container_poly &rcon, double &rx, double &ry, double &rz, double &rr, const double &ralfa, const double &rbeta, const double &rB, const double &riota)
{
	// [in]		con				the container with stored particles
	// [in]		x,y,z, r		point under examination
	// [in]     alfa, beta, B	feasibility constraints

	unsigned int i;
	int ijk_add, q_add, ijk, q;
	std::vector<int> cells;
	std::vector<int> neigh;
	voronoicell_neighbor c;

	// chci pridat bod, k tomu potrebuji nejake ID, jak ale zjistim, ktere ID je dostupne? co kdybych pouzil jako ID nejake
	//   cislo, ktere se bezne jako ID nepouziva (treba 0 nebo neco zaporneho)

	rcon.put(0, rx, ry, rz, rr);		// add particle under examination
	find_pos(ijk_add, q_add, 0, &rcon); // find position of added particle

//	rcon.compute_cell(c, ijk_add, q_add);
//	c.neighbors(neigh);

//  1) Feasibility for the neighborhood pr+sr:
//	cells.clear();
//
//	for (i = 0; i < neigh.size(); i++) {
//		find_pos(ijk, q, neigh[i], &rcon);
//		cells.push_back(ijk); cells.push_back(q); // add this neighbor to the vector
//	}
//	cells.push_back(ijk_add); cells.push_back(q_add);

//	if (feasibility(rcon, cells, ralfa, rbeta, rB)) { erase(ijk_add, q_add, &rcon); return 1; }

//  2) Feasibility for the whole container:
	if (overlap_f(rcon, riota)) {
		if (feasibility(rcon, ralfa, rbeta, rB)) { erase(ijk_add, q_add, &rcon); return 1; }
	}

	erase(ijk_add, q_add, &rcon);
	return 0;
}


// no_points returns number of points/generators in the subwindow determined by bounds lb, ub
int no_points(voro::container &rcon, double lb, double ub)
{
	// [in]		con		container with stored particles
	// [in]		lb, ub	bounds specifying the subwindow

	int i, j, id;
	int no = 0;
	double x, y, z;
	//int jo = 0;
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box
			id = rcon.id[j][i];
			x = rcon.p[j][3 * i];
			y = rcon.p[j][3 * i + 1];
			z = rcon.p[j][3 * i + 2];

			if (x > lb && x < ub && y > lb && y < ub && z > lb && z < ub) {
				// if the coordinates are in subwindow
				no++;
			} //else { jo++; }
		}
	}
	//std::cout << jo << " ";
	return no;
}

int no_points(voro::container_poly &rcon, double lb, double ub)
{
	// [in]		con		container with stored particles
	// [in]		lb, ub	bounds specifying the subwindow

	int i, j, id;
	int no = 0;
	double x, y, z;
	//int jo = 0;
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box
			id = rcon.id[j][i];
			x = rcon.p[j][4 * i];
			y = rcon.p[j][4 * i + 1];
			z = rcon.p[j][4 * i + 2];
			
			if (x > lb && x < ub && y > lb && y < ub && z > lb && z < ub) {
				// if the coordinates are in subwindow
				no++;
			} //else { jo++; }
		}
	}
	//std::cout << jo << " ";
	return no;
}

// only testing functions giving negative answer on the following question
// ??? muze opetovne pridana (tj. pridana po odebrani) castice zmenit box
int total_points(voro::container &rcon, const double &ralfa, const double &rbeta, const double &rB)
{
	int i, j, ijk, q, id;
	int no = 0;
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box
			id = rcon.id[j][0];
			find_pos(ijk, q, id, &rcon);	// find position of the particle
			try_delete(ijk, q, rcon, 1, ralfa, rbeta, rB);
			no++;
		}
	}
	return no;
}


// key function for estimating parameters theta and zet
void estim(double &th_estim, double &z_estim, voro::container &rcon, int N, double &ralfa_e, double &rbeta_e, double &rB_e, double lx, double ux, double ly, double uy, double lz, double uz)
{
	// [out]	th_estim					the estimate of parameter theta
	// [out]	z_estim						the estimate of parameter zet
	// [in]		con							the container with stored particles
	// [in]		N							the size of the sample for estimating integrals by Monte Carlo
	// [in]		alfa_e, beta_e, B_e			the estimates of hardcore parameters
	// [in]		lb, ub						lower and upper bound of subwindow

	int r, i;
	int na = 0;
	double eloc_e;
	bool add;
	double x,y,z;
	std::vector<int> va;
	std::vector<double> vb, vc;

	va.clear();
	vb.clear();
	vc.clear();

	// 0) ESTIM SUBWINDOW
	//estimates are done from a subwindow, i.e. removable points are computed on subwindow only
	// - bounds of subwindow are specified by doubles lb, ub (lb=0,ub=1 means to use the whole window)
	// - restriction on the subwindow helps to avoid boundary effects
	
	// 1A) no_removable returns number of removable points in container
	// - !! na rozdil od pridani bodu, ktere muze porusit horni mez, odebrani bodu dolni mez porusit nemuze
	//std::cout << "Number of removable: " << no_removable(rcon, 0, rbeta_e + 0.00001, rB_e, vb, 0, 1) << " from " << rcon.total_particles() << "\n";
	// - the question is whether use the correction rbeta_e + 0.00001 or not
	// - fc no_removable computes coefficients b(x)
	r = no_removable(rcon, ralfa_e, rbeta_e, rB_e, vb, lx, ux, ly, uy, lz, uz);
	if (r == 0) { std::cout << "ERROR: Estimates can not be computed (r)! \n"; return; } // ERROR when there are no removable points 
	else { std::cout << "Number of removable (subwindow " << lx << "-" << ux << ", " << ly << "-" << uy << ", " << lz << "-" << uz << "): " << r << "\n"; }
	 
	// alternative way how to compute number of removable points is usage of fc no_removable_2 (this fc was created for cross-checking)
	//r2 = no_removable_2(rcon, 0, rbeta_e + 0.00001, rB_e);
	//std::cout << "Number of removable (2): " << r2 << " from " << rcon.total_particles() << "\n";
	
	// 1B) feasibility under estimates 
	// - it is important to verify if estimated values of hardcore parameters give the feasibility of simulated configuration (if not they have to be corrected somehow)
	// - ?? but how to correct ??
	if (feasibility(rcon, ralfa_e, rbeta_e, rB_e) == 1) { std::cout << "Feasible under estimates. \n"; }
	else { std::cout << "Not feasible under estimates! \n"; }

	// 1C) vypocet integralu - koeficientu a(x), c(x) (PLL_coeff)
	// 1C.I) using the sample from the uniform distribution
	for (i = 0; i < N; i++) {
		// generating of the sample (coordinates from (lb,ub)) from uniform distribution:
		x = lx + (ux-lx)*uniform(0, 1); y = ly + (uy-ly)*uniform(0, 1); z = lz + (uz-lz)*uniform(0, 1);  

		add = addable(rcon, x, y, z, ralfa_e, rbeta_e, rB_e); // try if the particle is addable or not (returns 0 if not, 1 if yes)
		//add = addable(rcon, x, y, z, 0, ralfa_e, rB_e); // it is important to have feasible starting configuration

		if(add==true){ va.push_back(1); } else { va.push_back(0); } // save this information

		if (add == 0) {
			// loc energy is infite = point is not addable
			vc.push_back(0); // save arbitrary value - in this case the value is not important
		}
		else {
			eloc_e = try_add(0, x, y, z, rcon, 1, ralfa_e, rbeta_e, rB_e); // compute the adding energy iff addable (for theta value equal 1)
			vc.push_back(-log(eloc_e)); // save the correct value
		}
		// !!!!!! pozor na pripustnost s beta_e

		// fc try_add returns zero iff the point is not addable
		//    therefore we do not necessarily need to use fc addable, but on the other hand it is
		//	  a little bit unreliable to test if double == 0
	}
	for (i = 0; i < N; i++) { na = na + va[i]; }
	if (na == 0) { std::cout << "ERROR: Estimates can not be computed (a)! \n"; return; } // ERROR when there are no addable points (i.e. integral can not be computed)
	std::cout << "Number of addable: " << na << " from " << N << "\n";

	// 1C.II) using the grid
/*	int s = 20;
	N = s*s*s;
	int j, k, l;
	int p = 1 / (2 * s);
	for (j = 0; j < s; j++) {
		x = p + j / s;
		for (k = 0; k < s; k++) {
			y = p + k / s;
			for (l = 0; l < s; l++) {
				z = p + l / s;
				add = addable(rcon, x, y, z, ralfa_e, rbeta_e, rB_e); // try if the particle is addable or not (returns 0 if not, 1 if yes)
																	 
				if (add == true) { va.push_back(1); }
				else { va.push_back(0); } // save this information

				if (add == 0) {
					// loc energy is infite = point is not addable
					vc.push_back(0); // save arbitrary value - in this case the value is not important
				}
				else {
					eloc_e = try_add(0, x, y, z, rcon, 1, ralfa_e, rbeta_e, rB_e); // compute the adding energy iff addable (for theta value equal 1)
					vc.push_back(-log(eloc_e)); // save the correct value
				}
			}
		}
	}
	for (i = 0; i < N; i++) { na = na + va[i]; }
	if (na == 0) { std::cout << "ERROR: Estimates can not be computed (a-g)! \n"; return; } // ERROR when there are no addable points (i.e. integral can not be computed)
	std::cout << "Number of addable (grid): " << na << " from " << N << "\n";
*/
	// vypocet sumy - koeficientu b(x) (PLL_const)
	//  ... uvnitr fce no_removable
	std::cout << "\n";
	
	// 2) Takacz-Fiksel estimate can be obtained for different choices of testing functions h1 and h2
	// 2A)------------------------------------------------------ Quick estimate --------------------------------------------------------------------------
	// - first testing fc is identically equal to 1

	//bisection(va, vb, vc, k, 0);							// 1) bisection method (under construction)
	//secant(va, vb, vc, k, 0);								// 2) secant method (under construction)
	th_estim = NR(va, vb, vc, r, N);						// 3) Newton-Raphson method

	z = 0;
	for (i = 0; i < N; i++) {
		z = z + exp(-th_estim*vc[i])*va[i];
	}
	z_estim = (N*r) / z;
	

	std::cout << "Q: Estimates of (theta,z): " << th_estim << " " << z_estim << "\n";


	// 2B)---------------------------------------------------- Max. pseudolikelihood estimate ------------------------------------------------------------
	// 2B.I) classic approach searching for roots of derivation equation
	
	//bisection(va, vb, vc, k, 0);							// 1) bisection method
	//secant(va, vb, vc, k, 0);								// 2) secant method
	th_estim = NR(va, vb, vc, r, 0);						// 3) Newton-Raphson method

	z = 0;
	for (i = 0; i < N; i++) {
		z = z + exp(-th_estim*vc[i])*va[i];
	}
	z_estim = (N*r) / z;
	// zet_estim(N, th_estim, va, vc); // estimation of intesnsity parameter --> testing fc is identically equal 1

	std::cout << "MPLE-c: Estimates of (theta,z): " << th_estim << " " << z_estim << "\n";

	// 2B.II) approach using LBFGS optimization algorithm directly to pseudolikelihood fc insted of its derivatives
	// - uses the library LBFGS++ described here: https://yixuan.cos.name/LBFGSpp/	
	/*
	const int n = 10;
	// Set up parameters
	LBFGSpp::LBFGSParam<double> param;
	param.epsilon = 1e-6;
	param.max_iterations = 100;

	// Create solver and function object
	LBFGSpp::LBFGSSolver<double> solver(param);
	PLL fun(n, va, vb, vc, r);

	// Initial guess
	Eigen::VectorXd ex = Eigen::VectorXd::Zero(n);
	// x will be overwritten to be the best point found
	double fx;
	int niter = solver.minimize(fun, ex, fx);

	std::cout << niter << " iterations" << std::endl;
	std::cout << "x = \n" << ex.transpose() << std::endl;
	std::cout << "f(x) = " << fx << std::endl;
	
	// ?? x contains estimates of theta and zet --> how?
	//th_estim = x[0]; z_estim = x[1];
	//std::cout << "MPLE-LBFGS: Estimates of (theta,z): " <<< th_estim << " " << z_estim << "\n";
	//std::cout << "            Number of iterations was " << niter << " and optimal value reached was " << fx << "\n";
	
*/

	// 3) using an initial guess on parameters (Takacz-Fiksel estimate can be used) the maximum likelihood procedure can be started
	// - improves the estimate values from the initial guess
	// - function MLE under construction below

}




void estim(double &th_estim, double &z_estim, voro::container_poly &rcon, int N, double &ralfa_e, double &rbeta_e, double &rB_e, double &riota_e, double lx, double ux, double ly, double uy, double lz, double uz)
{
	// [out]	th_estim					the estimate of parameter theta
	// [out]	z_estim						the estimate of parameter zet
	// [in]		con							the container with stored particles containing the particles radii
	// [in]		N							the size of the sample for estimating integrals by Monte Carlo
	// [in]		alfa_e, beta_e, B_e			the estimates of hardcore parameters
	// [in]		lb, ub						lower and upper bound of subwindow (different bounds in every coordinate)

	int r, i;
	int na = 0;
	double eloc_e;
	bool add;
	double x, y, z, rr;
	std::vector<int> va;
	std::vector<double> vb, vc;

	va.clear();
	vb.clear();
	vc.clear();

	// 0) ESTIM SUBWINDOW
	//estimates are done from a subwindow, i.e. removable points are computed on subwindow only
	// - bounds of subwindow are specified by doubles lb, ub (lb=0,ub=1 means to use the whole window)
	// - restriction on the subwindow helps to avoid boundary effects

	// 1A) no_removable returns number of removable points in container
	// - !! na rozdil od pridani bodu, ktere muze porusit horni mez, odebrani bodu dolni mez porusit nemuze
	//std::cout << "Number of removable: " << no_removable(rcon, 0, rbeta_e + 0.00001, rB_e, vb, 0, 1) << " from " << rcon.total_particles() << "\n";
	// - the question is whether use the correction rbeta_e + 0.00001 or not
	// - fc no_removable computes coefficients b(x)
	r = no_removable(rcon, 0, rbeta_e, rB_e, riota_e, vb, lx, ux, ly, uy, lz, uz);
	if (r == 0) { std::cout << "ERROR: Estimates can not be computed (r)! \n"; return; }
	else { std::cout << "Number of removable (subwindow): " << r << "\n"; }

	// 1B) vypocet integralu - koeficientu a(x), c(x) (PLL_coeff)
	for (i = 0; i < N; i++) {
		x = lx + (ux - lx)*uniform(0, 1); y = ly + (uy - ly)*uniform(0, 1); z = lz + (uz - lz)*uniform(0, 1);  
		//rr = triangle(0.005, 0.03, 0.02); // generating of the sample - coordinates from (lb,ub) and radius from ... 
		//rr = 0.005 + 0.025*uniform(0, 1); 
		rr = 0.05*uniform(0, 1);
		//rr = gamma(3, 0.5);
		// or using the average radius in container:
		//rr = ave_rad(conp);
		//rr = 100*uniform(0,1);
		// (the same sampler for radius as in try_add fc is used)

		add = addable(rcon, x, y, z, rr, ralfa_e, rbeta_e, rB_e, riota_e); // try if the particle is addable or not
		va.push_back(add); // save this information

		if (add == 0) {
			// loc energy is infite = point is not addable
			vc.push_back(0); // save arbitrary value - in this case the value is not important
		}
		else {
			eloc_e = try_add(0, x, y, z, rr, rcon, 1, ralfa_e, rbeta_e, rB_e, riota_e); // compute the adding energy iff addable
			vc.push_back(-log(eloc_e)); // save the correct value
		}
	}

	for (i = 0; i < N; i++) { na = na + va[i]; }
	if (na == 0) { std::cout << "ERROR: Estimates can not be computed (a)! \n"; return; }
	std::cout << "Number of addable: " << na << " from " << N << "\n";

	std::cout << "\n";

	// 2) Takacz-Fiksel estimate can be obtained for different choices of testing functions h1 and h2
	// 2A)------------------------------------------------------ Quick estimate --------------------------------------------------------------------------
	// - first testing fc is identically equal to 1
	//bisection(va, vb, vc, k, 0);
	//secant(va, vb, vc, k, 0);
	th_estim = NR(va, vb, vc, r, N);		// 3) Newton-Raphson method

	z = 0;
	for (i = 0; i < N; i++) {
		z = z + exp(-th_estim*vc[i])*va[i];
	}
	z_estim = (N*r) / z;
	// zet_estim(N, th_estim, va, vc); // estimation of intesnsity parameter --> testing fc is identically equal 1

	std::cout << "Q: Estimates of (theta,z): " << th_estim << " " << z_estim << "\n";


	// 2B)---------------------------------------------------- Max. pseudolikelihood estimate ------------------------------------------------------------
	// 2B.I) classic approach searching for roots of derivation equation

	//bisection(va, vb, vc, k, 0);							// 1) bisection method
	//secant(va, vb, vc, k, 0);								// 2) secant method
	th_estim = NR(va, vb, vc, r, 0);						// 3) Newton-Raphson method

	z = 0;
	for (i = 0; i < N; i++) {
		z = z + exp(-th_estim*vc[i])*va[i];
	}
	z_estim = (N*r) / z;
	// zet_estim(N, th_estim, va, vc); // estimation of intesnsity parameter --> testing fc is identically equal 1

	std::cout << "MPLE-c: Estimates of (theta,z): " << th_estim << " " << z_estim << "\n";

	// 2B.II) approach using LBFGS optimization algorithm directly to pseudolikelihood fc insted of its derivatives
	// - uses the library LBFGS++ described here: https://yixuan.cos.name/LBFGSpp/	
	/*
	//	0) LBFGS method
	const int n = 10;
	// Set up parameters
	LBFGSpp::LBFGSParam<double> param;
	param.epsilon = 1e-6;
	param.max_iterations = 100;

	// Create solver and function object
	LBFGSpp::LBFGSSolver<double> solver(param);
	PLL fun(n, va, vb, vc, r);

	// Initial guess
	Eigen::VectorXd ex = Eigen::VectorXd::Zero(n);
	// x will be overwritten to be the best point found
	double fx;
	int niter = solver.minimize(fun, ex, fx);

	std::cout << niter << " iterations" << std::endl;
	std::cout << " x = \n" << ex.transpose() << std::endl;
	std::cout << " f(x) = " << fx << std::endl;

	// ?? x contains estimates of theta and zet --> how?
	//th_estim = x[0]; z_estim = x[1];
	//std::cout << "MPLE-LBFGS: Estimates of (theta,z): " << th_estim << " " << z_estim << "\n";
	//std::cout << "            Number of iterations was " << niter << " and optimal value reached was " << fx << "\n";
	*/


	// 3) using an initial guess on parameters (Takacz-Fiksel estimate can be used) the maximum likelihood procedure can be started
	// - improves the estimate values from the initial guess
	// - function MLE under construction below

	
}



double zet_estim(int rN, double &rth, std::vector<bool> &ra, std::vector<double> &rc)
{
	// [in]		N				the size of the sample for estimating integrals by Monte Carlo
	// [in]		th				the estimate of parameter theta
	// [in]		a				the vector of 0,1 saying if the i-th point is addable or not
	// [in]		c				the vector of local energies for sampled points
	int i;
	double z = 0;

	// approximation = plain Monte Carlo method
	for (i = 0; i < rN; i++) {

		z = z + exp(-rth*rc[i])*ra[i];

	}

	return z;
}



// fc computing the constant part of the pseudolikelihood contrast function
double PLL_const(voro::container &rcon, double &ralfa_e, double &rbeta_e, double &rB_e)
{
	// [in]		con				the container with stored particles
	// [in]		alfa_e			the estimate of hardcore parameter alfa
	// [in]		beta_e			the estimate of hardcore parameter beta
	// [in]		B_e				the estimate of hardcore parameter B

	double k = 0;
	double eloc_e;
	int i, j;
	int ijk,q;

	//	c_loop_all clo(rcon);          // loop pres vsechny castice - nelze pouzit, prtz fc try_delete meni strukturu container
	// alternativa: loop pres vsechna ID

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box - prtz se castice smaze a pak zase prida,
										   // je potreba co[j]-krat vzit prvni castici v bodu, ta se vzdy presune na konec
			 // musi byt ijk == j, q == 0 
			/*id = rcon.id[j][0];
			find_pos(ijk, q, id, &rcon);
			std::cout << ijk << " " << q << " ID: " << id << " " << rcon.id[ijk][q] << "\n";
			*/
			ijk = j; q = 0; // je potreba aktualizovat kazdy cyklus, nebot fce try_delete to zmeni
			// I can try if the point uder consideration is removable, if yes I can delete it, call fc try_add which computes
			//		exp of local energy and after that I can add the point back to the container
			// or
			// Instead finding out if it is removable I can call fc try_del, it will verify removeability (if not it returns 0)
			//		and moreover it computes the value of local energy but in the opposite direction, thus the local energy
			//		will be equal to -ln(try_del) instead of ln(try_add)

			// !! POZOR: try_delete zmeni ijk,q
			eloc_e = try_delete(ijk, q, rcon, 1, ralfa_e, rbeta_e, rB_e);

			if (eloc_e != 0) {
				// s um only over removable points - if the point cant be removed, fc try_delete returns 0
				//	(otherwise it returns positive number - exp is positive)
				k = k + log(eloc_e); 
			}
		}
	}

	//	} while (clo.inc());    // inc() - finds the next particle to test.

	return k;
}


// fc computing the coefficients of the pseudolikelihood contrast function
//   note: size of the sample (N) should not be fixed, but it should correspond to accuracy 
//		(N should be so big so that the variance was smaller than a given value) - fc then returns N
long PLL_coeff(voro::container &rcon, double &ralfa_e, double &rbeta_e, double &rB_e, int rN, std::vector<bool> &ra, std::vector<double> &rc)
{
	// [in]		con				the container with stored particles
	// [in]		alfa_e			the estimate of hardcore parameter alfa
	// [in]		beta_e			the estimate of hardcore parameter beta
	// [in]		N				the size of the sample for estimating integrals by Monte Carlo
	// [out]	a				the vector of 0,1 saying if the i-th point is addable or not
	// [out]	c				the vector of local energies for sampled points

	int i;
	double x, y, z;
	double eloc_e;
	bool add;

	ra.clear();
	rc.clear();

	for (i = 0; i < rN; i++) {
		x = uniform(0, 1); y = uniform(0, 1); z = uniform(0, 1);  // generating of the sample

		add = addable(rcon, x, y, z, ralfa_e, rbeta_e, rB_e); // try if the particle is addable or not
		ra.push_back(add); // save this information

		if (add == 0) {
			// loc energy is infite = point is not addable
			rc.push_back(0); // save arbitrary value - in this case the value is not important
		}
		else {
			eloc_e = try_add(0, x, y, z, rcon, 1, ralfa_e, rbeta_e, rB_e); // compute the adding energy iff addable
			rc.push_back(-log(eloc_e)); // save the correct value
		}

		// fc try_add returns zero iff the point is not addable
		//    therefore we do not necessarily need to use fc addable, but on the other hand it is
		//	  a little bit unreliable to test if double == 0

	}    
	return 1;
}




void raw_resid(std::vector<double> &raw_res, voro::container &rcon, int &gsi, int rN, double &ralfa_e, double &rbeta_e, double &rB_e, double &th_estim, double &z_estim)
{
	// [out]	raw_res			vector of raw residulas (residaul computed for every cube)
	// [in]		con				the container with stored particles
	// [in]		N				the size of the sample for estimating integrals by Monte Carlo
	// [in]		gsi				the grid size (number of cubes in each direction)
	// [in]		alfa_e			the estimate of hardcore parameter alfa
	// [in]		beta_e			the estimate of hardcore parameter beta
	// [in]		B_e				the estimate of hardcore parameter B
	// [in]		th_estim		the estimate of parameter theta
	// [in]		z_estim			the estimate of parameter zet

	int i, j, k, l;
	int x, y, z, del;
	double xx, yy, zz;
	double eloc_e = 0;

	raw_res.clear();				//vycistit
	raw_res.reserve(gsi*gsi*gsi);	//naalokovat pamet (vime kolik budeme potrebovat)
	// computing # REMOVABLE POINTS
	// set the null vector of length gsi^3:
	std::vector<int> nor(gsi*gsi*gsi); 
	// loop over all generators -- mod its coordinates -- obtain j,k,l values -- inc appropriate value in the vector
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			del = rcon.id[j][0]; // fc removable changes structure
			// bez teto podminky if, by se rozpocitaly vsechny body --> kontrola: suma vektoru raw_res by pak mela dat pocet castic v containeru (viz point_density)
			if (removable(rcon, del, ralfa_e, rbeta_e, rB_e)) {

				x = (int)(rcon.p[j][0] * gsi);      // coordinates --> i,j,k values
				y = (int)(rcon.p[j][1] * gsi);		// (int) : 0.897*10 -> 8 , 0.210*50 -> 10  ...  dava poradi mnoziny do ktere spadne, cislovano od 0
				z = (int)(rcon.p[j][2] * gsi);

				(nor[x*gsi*gsi + y*gsi + z])++;
			}
		}
	}

	// computing RAW RESIDUALS
	for (j = 0; j < gsi; j++) {
		for (k = 0; k < gsi; k++) {
			for (l = 0; l < gsi; l++) {
				
				for (i = 0; i < rN; i++) { // gsi*gsi*gsi*N-krat se vola try_add (napr. gsi = 10, N = 1000 ----> 1 000 000 !!!!!!!!! to je dost !!!!!!!!
					//xx = uniform(j / gsi, (j + 1) / gsi); yy = uniform(k / gsi, (k + 1) / gsi); zz = uniform(l / gsi, (l + 1) / gsi);  // generating of the sample
					xx = j / gsi + uniform(0, 1) / gsi; yy = k / gsi + uniform(0, 1) / gsi; zz = l / gsi + uniform(0, 1) / gsi;  // generating of the sample

					eloc_e = eloc_e + try_add(0, xx, yy, zz, rcon, th_estim, ralfa_e, rbeta_e, rB_e); // sum up the exp of adding energy (0 iff not addable)

				}
				eloc_e = eloc_e / rN;
				raw_res.push_back(nor[j*gsi*gsi + k*gsi + l] - z_estim * eloc_e); // save residual
																				// nor is vector of removable points per cubes
			}
		}
	}

	// other types of residuals --> another functions
}

void raw_resid(std::vector<double> &raw_res, voro::container_poly &rcon, int &gsi, int rN, double &ralfa_e, double &rbeta_e, double &rB_e, double &riota_e, double &th_estim, double &z_estim)
{
	// [out]	raw_res			vector of raw residulas (residaul computed for every cube)
	// [in]		con				the container with stored particles
	// [in]		N				the size of the sample for estimating integrals by Monte Carlo
	// [in]		gsi				the grid size (number of cubes in each direction)
	// [in]		alfa_e			the estimate of hardcore parameter alfa
	// [in]		beta_e			the estimate of hardcore parameter beta
	// [in]		B_e				the estimate of hardcore parameter B
	// [in]		th_estim		the estimate of parameter theta
	// [in]		z_estim			the estimate of parameter zet

	int i, j, k, l;
	int x, y, z, del;
	double xx, yy, zz, rr;
	double eloc_e = 0;

	raw_res.clear();				//vycistit
	raw_res.reserve(gsi*gsi*gsi);	//naalokovat pamet (vime kolik budeme potrebovat)
	// computing # REMOVABLE POINTS
	// set the null vector of length gsi^3:
	std::vector<int> nor(gsi*gsi*gsi);
	// loop over all generators -- mod its coordinates -- obtain j,k,l values -- inc appropriate value in the vector
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			del = rcon.id[j][0]; // fc removable changes structure
								 // bez teto podminky if, by se rozpocitaly vsechny body --> kontrola: suma vektoru raw_res by pak mela dat pocet castic v containeru (viz point_density)
			if (removable(rcon, del, ralfa_e, rbeta_e, rB_e)) {

				x = (int)(rcon.p[j][0] * gsi);      // coordinates --> i,j,k values
				y = (int)(rcon.p[j][1] * gsi);		// (int) : 0.897*10 -> 8 , 0.210*50 -> 10  ...  dava poradi mnoziny do ktere spadne, cislovano od 0
				z = (int)(rcon.p[j][2] * gsi);

				(nor[x*gsi*gsi + y*gsi + z])++;
			}
		}
	}

	// computing RAW RESIDUALS
	for (j = 0; j < gsi; j++) {
		for (k = 0; k < gsi; k++) {
			for (l = 0; l < gsi; l++) {

				for (i = 0; i < rN; i++) { // gsi*gsi*gsi*N-krat se vola try_add (napr. gsi = 10, N = 1000 ----> 1 000 000 !!!!!!!!! to je dost !!!!!!!!
										   //xx = uniform(j / gsi, (j + 1) / gsi); yy = uniform(k / gsi, (k + 1) / gsi); zz = uniform(l / gsi, (l + 1) / gsi);  // generating of the sample
					xx = j / gsi + uniform(0, 1) / gsi; yy = k / gsi + uniform(0, 1) / gsi; zz = l / gsi + uniform(0, 1) / gsi;  // generating of the sample
					rr = triangle(0.005, 0.03, 0.02);

					eloc_e = eloc_e + try_add(0, xx, yy, zz, rr, rcon, th_estim, ralfa_e, rbeta_e, rB_e, riota_e); // sum up the exp of adding energy (0 iff not addable)

				}
				eloc_e = eloc_e / rN;
				raw_res.push_back(nor[j*gsi*gsi + k*gsi + l] - z_estim * eloc_e); // save residual
																				  // nor is vector of removable points per cubes
			}
		}
	}

	// other types of residuals --> another functions
}


void inverse_resid(std::vector<double> &inv_res, voro::container &rcon, int &gsi, int rN, double &ralfa_e, double &rbeta_e, double &rB_e, double &th_estim, double &z_estim)
{
	// [out]	inv_res			vector of inverse residulas (residaul computed for every cube)
	// [in]		con				the container with stored particles
	// [in]		N				the size of the sample for estimating integrals by Monte Carlo
	// [in]		gsi				the grid size (number of cubes in each direction)
	// [in]		alfa_e			the estimate of hardcore parameter alfa
	// [in]		beta_e			the estimate of hardcore parameter beta
	// [in]		B_e				the estimate of hardcore parameter B
	// [in]		th_estim		the estimate of parameter theta
	// [in]		z_estim			the estimate of parameter zet

	int i, j, k, l;
	int x, y, z, del;
	double xx, yy, zz;
	double eloc_e;
	int nul = 0;

	inv_res.clear();				//vycistit
	inv_res.reserve(gsi*gsi*gsi);	//naalokovat pamet (vime kolik budeme potrebovat)
	// computing gain of REMOVABLE POINTS
	// set the null vector of length gsi^3:
	std::vector<double> nor(gsi*gsi*gsi);
	// loop over all generators -- mod its coordinates -- obtain j,k,l values -- inc appropriate value in the vector
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			del = rcon.id[j][0]; // fc removable changes structure
			if (removable(rcon, del, ralfa_e, rbeta_e, rB_e)) {

				// potrebuju spocitat exp -energie pridani odstranitelne castice == exp energie odebrani teto castice == exp(-log(try_delete))
				// try_delete == exp -energie odebrani castice ; -energie pridani == energie odebrani
				eloc_e = try_delete(j, nul, rcon, 1, ralfa_e, rbeta_e, rB_e); 
				// podminka if removable je zde proto abych nemusel zkouset if eloc_e > 0 (double > 0)
				eloc_e = exp(-log(eloc_e));

				x = (int)(rcon.p[j][0] * gsi);      // coordinates --> i,j,k values
				y = (int)(rcon.p[j][1] * gsi);		// (int) : 0.897*10 -> 8 , 0.210*50 -> 10  ...  dava poradi mnoziny do ktere spadne, cislovano od 0
				z = (int)(rcon.p[j][2] * gsi);

				nor[x*gsi*gsi + y*gsi + z] = nor[x*gsi*gsi + y*gsi + z] + 1/eloc_e;
			}
		}
	}

	eloc_e = 0;
	// computing INVERSE RESIDUALS
	for (j = 0; j < gsi; j++) {
		for (k = 0; k < gsi; k++) {
			for (l = 0; l < gsi; l++) {

				for (i = 0; i < rN; i++) {
					//xx = uniform(j / gsi, (j + 1) / gsi); yy = uniform(k / gsi, (k + 1) / gsi); zz = uniform(l / gsi, (l + 1) / gsi);  // generating of the sample
					xx = j/gsi + uniform(0, 1)/gsi; yy = k/gsi + uniform(0, 1)/gsi; zz = l/gsi + uniform(0, 1)/gsi;  // generating of the sample

					//eloc_e = try_add(0, xx, yy, zz, rcon, th_estim, ralfa_e, rbeta_e, rB_e); // sum up the exp of adding energy (0 iff not addable)
					if (addable(rcon, xx, yy, zz, ralfa_e, rbeta_e, rB_e)) { nul++; }
				}
				nul = nul / rN;
				inv_res.push_back(nor[j*gsi*gsi + k*gsi + l] - z_estim * nul); // save residual
																				  // nor is vector of gains of removable points per cubes
			}
		}
	}

	// other types of residuals --> another functions
}


void pearson_resid(std::vector<double> &p_res, voro::container &rcon, int &gsi, int rN, double &ralfa_e, double &rbeta_e, double &rB_e, double &th_estim, double &z_estim)
{
	// [out]	p_res			vector of pearson residulas (residaul computed for every cube)
	// [in]		con				the container with stored particles
	// [in]		N				the size of the sample for estimating integrals by Monte Carlo
	// [in]		gsi				the grid size (number of cubes in each direction)
	// [in]		alfa_e			the estimate of hardcore parameter alfa
	// [in]		beta_e			the estimate of hardcore parameter beta
	// [in]		B_e				the estimate of hardcore parameter B
	// [in]		th_estim		the estimate of parameter theta
	// [in]		z_estim			the estimate of parameter zet

	int i, j, k, l;
	int x, y, z, del;
	double xx, yy, zz;
	double eloc_e;
	int nul = 0;

	p_res.clear();				//vycistit
	p_res.reserve(gsi*gsi*gsi);	//naalokovat pamet (vime kolik budeme potrebovat)
	// computing gain of REMOVABLE POINTS
	// set the null vector of length gsi^3:
	std::vector<double> nor(gsi*gsi*gsi);
	// loop over all generators -- mod its coordinates -- obtain j,k,l values -- inc appropriate value in the vector
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			del = rcon.id[j][0]; // fc removable changes structure
			if (removable(rcon, del, ralfa_e, rbeta_e, rB_e)) {

				// potrebuju spocitat exp -energie pridani odstranitelne castice == exp energie odebrani teto castice == exp(-log(try_delete))
				// try_delete == exp -energie odebrani castice ; -energie pridani == energie odebrani
				eloc_e = try_delete(j, nul, rcon, 1, ralfa_e, rbeta_e, rB_e);
				// podminka if removable je zde proto abych nemusel zkouset if eloc_e > 0 (double > 0)
				eloc_e = exp(-log(eloc_e));

				x = (int)(rcon.p[j][0] * gsi);      // coordinates --> i,j,k values
				y = (int)(rcon.p[j][1] * gsi);		// (int) : 0.897*10 -> 8 , 0.210*50 -> 10  ...  dava poradi mnoziny do ktere spadne, cislovano od 0
				z = (int)(rcon.p[j][2] * gsi);

				nor[x*gsi*gsi + y*gsi + z] = nor[x*gsi*gsi + y*gsi + z] + 1 / sqrt(eloc_e);
			}
		}
	}

	eloc_e = 0;
	// computing PEARSON RESIDUALS
	for (j = 0; j < gsi; j++) {
		for (k = 0; k < gsi; k++) {
			for (l = 0; l < gsi; l++) {

				for (i = 0; i < rN; i++) {
					//xx = uniform(j / gsi, (j + 1) / gsi); yy = uniform(k / gsi, (k + 1) / gsi); zz = uniform(l / gsi, (l + 1) / gsi);  // generating of the sample
					xx = j / gsi + uniform(0, 1) / gsi; yy = k / gsi + uniform(0, 1) / gsi; zz = l / gsi + uniform(0, 1) / gsi;  // generating of the sample

					eloc_e = eloc_e + sqrt(try_add(0, xx, yy, zz, rcon, th_estim, ralfa_e, rbeta_e, rB_e)); // sum up the exp of sqrt of adding energy (0 iff not addable)

				}
				eloc_e = eloc_e / rN;
				p_res.push_back(nor[j*gsi*gsi + k*gsi + l] - z_estim * eloc_e); // save residual
																				  // nor is vector of gains of removable points per cubes
			}
		}
	}

	// other types of residuals --> another functions
}


void smooth_res_field(voro::container &rcon, double &ralfa_e, double &rbeta_e, double &rB_e, double &th_estim, double &z_estim)
{

}


//fc MLE computes max likelihood estimator for parameter theta; before computations its necessary to specify:
//			- theta; alfa, beta, B, z ... initial value to be improved; fixed values of parameters estimated before
//			- n ... number of simulations
//			- it ... number of iterations per simulation
//			- lb, ub ... bounds of used subwindow
//			- pcon ... data to be imported
//			- critical value, which when overcomed, it is necessary to do new simulations
//			- stopping condition (to specify when the estimate is good enough) - use rather accuracy condition or rather specify a given number of iterations?
void MLE(double th_estim, double alfa_e, double beta_e, double B_e, double zet_e)
{
	// [in,out]	th_estim	the initial value is iteratively improved
	// [in]		alfa_e,...	estimated parameters considered to be fixed

	int i, j;
	long no;
	int nx, ny, nz;
	double theta_e = th_estim;
	double theta_e_prev;

	// set up
	int n = 10;	// the number of simulations
	std::vector<double> t(n); // the vector for storing the values of sufficient statistic t(Y_0), ..., t(Y_n-1)
	int it = 200000; // number of iterations per simulation
	long n_add = 0; // citace poctu pridanych/smazanych/posunutych castic (celkovy pocet castic lze odtud pak dopocitat)
	long n_del = 0;						// prebytecne, nechci tu zaznam run 
	long n_mov = 0;
	//long n_part;
	double lb = 0.2; // bounds of subwindow
	double ub = 0.8;
	std::vector<int> fid;
	fid.clear();
	double sigma = pow(0.015, 2);
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;
	double p;
	double e = 0.01; // stopping condition on accuracy
	int s = 15; // number of steps specification
	int counter = 0;
	double d = 0.2; // critical value for new simulations

	// import of data
	pre_container pcon(0, 1, 0, 1, 0, 1, true, true, true);  // true = periodic in given coordinate
	pcon.import("C://Users/Vigo/Dokumenty/Visual Studio 2015/Projects/Gibbs-Voronoi/Gibbs-Voronoi/SimulationsPRPAP/study4/src1.0/datacon.txt");
	pcon.guess_optimal(nx, ny, nz);  // guess
	container con(0, 1, 0, 1, 0, 1, nx, ny, nz, true, true, true, 8);

	double tt = 0; // energie puvodnich dat -- nebylo by lepsi predat fci parametrem???

	// 1st simulations
	for (i = 0; i < n; i++) {
		con.clear();
		pcon.setup(con);  // the setting of container and import of particles
		no = con.total_particles();
		fid.clear();

		for (j = 0; j < it; j++) {
			//std::cout << " STEP " << i++ << '\n';
			bdma_step(no, fid, con, sigma, th_estim, zet_e, alfa_e, beta_e, B_e, n_add, n_del, n_mov);
		}

		t[i] = con_energy(con, lb, ub); // energy of the container - pouzit podokno???
	}

	// zastitit nejakym while cyklem ... dokud zmena theta_e nebude dostatecne mala (ve vicero po sobe jdoucich krocich)
	do
	{
		counter++;
		// update
		for (i = 0; i < n; i++) {
			sum1 = sum1 + (exp(-theta_e*t[i] + th_estim*t[i]));
			sum2 = sum2 + t[i] * (exp(-theta_e*t[i] + th_estim*t[i]));
			sum3 = sum3 + pow(t[i], 2)*(exp(-theta_e*t[i] + th_estim*t[i]));
		}
		p = sum2 / sum1;

		theta_e_prev = theta_e;
		theta_e = theta_e + (-tt + p) / (sum3 - pow(p, 2));

		// new simulations
		if (abs(theta_e - th_estim) > d) {
			/*do a new simulations*/ 
			th_estim = theta_e; // ??? mela by se aktualizovat inicializacni hodnota theta_0 (tj. melo by dojit ke zmene v sumach) ???
								
			for (i = 0; i < n; i++) {
				con.clear();
				pcon.setup(con);  // the setting of container and import of particles
				no = con.total_particles();
				fid.clear();

				for (j = 0; j < it; j++) {
					//std::cout << " STEP " << i++ << '\n';
					bdma_step(no, fid, con, sigma, th_estim, zet_e, alfa_e, beta_e, B_e, n_add, n_del, n_mov);
				}

				t[i] = con_energy(con, lb, ub); // energy of the container - pouzit podokno???
			}
		}
	//} while (abs(theta_e - theta_e_prev) > e);  // accuracy condition
	} while (counter < s);	// specification of the number of steps

}

// con_energy returns the average value of the part of energy which does not depend on estimated parameter
// ??? lepe prepsat nasledujici fci jako metodu tridy ???
double con_energy(voro::container &rcon, double lb, double ub)
{
	// [in]		con		container with stored particles
	// [in]		lb,ub	bounds of the subwindow

	int i, j, id, ijk, q;
	unsigned int k, fng;
	double t = 0;
	double x, y, z, xn, yn, zn;

	int cit = 0;
	double enr;

	voronoicell_neighbor c,d;  // bunka
	std::vector<int> neigh;  // vektor pro prirazeni ID sousedu
	std::vector<int> vert;

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box
			x = rcon.p[j][3 * i];
			y = rcon.p[j][3 * i + 1];
			z = rcon.p[j][3 * i + 2];
			if (x > lb && x < ub && y > lb && y < ub && z > lb && z < ub) {
				// if the coordinates are in subwindow
				id = rcon.id[j][i];
				rcon.compute_cell(c, j, i);
				c.neighbors(neigh);  // Computes a vector list of neighbors of the current particle under consideration
				c.face_vertices(vert);
				fng = 0;
									 // loop over the neighbors
				for (k = 0; k < neigh.size(); k++) {

					// Skip if the neighbor information is smaller than this particle's ID, to avoid double counting. This also
					// removes faces that touch the walls, since the neighbor information is set to negative numbers for these cases.
					if (neigh[k] > id) {

						find_pos(ijk, q, neigh[k], &rcon);  // find position of the i-th neighbor

						face_dist(fng, vert, x, y, z, xn, yn, zn, c);		// {1PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
						xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;	// {1PER} urcim skutecne souradnice tohoto souseda za stenou

						rcon.compute_cell(d, ijk, q);	// compute the cell of particle
						enr = V2(c, d, x, y, z, xn, yn, zn);
						t = t + enr;				// increase energy
						if (enr > 0) { cit++; }		// increase number of pairs creating the energy (only if energy is > 0)

					}
					fng = fng + vert[fng] + 1;
				}
			}
			
		}
	}

	//return t;
	return t / cit;
}

double con_energy(voro::container_poly &rcon, double lb, double ub)
{
	// [in]		con		container with stored particles
	// [in]		lb,ub	bounds of the subwindow

	int i, j, id, ijk, q;
	unsigned int k, fng;
	double t = 0;
	double x, y, z, xn, yn, zn;

	int cit = 0;
	double enr;

	voronoicell_neighbor c, d;  // bunka
	std::vector<int> neigh;  // vektor pro prirazeni ID sousedu
	std::vector<int> vert;

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box
			x = rcon.p[j][4 * i];
			y = rcon.p[j][4 * i + 1];
			z = rcon.p[j][4 * i + 2];
			if (x > lb && x < ub && y > lb && y < ub && z > lb && z < ub) {
				// if the coordinates are in subwindow
				id = rcon.id[j][i];
				rcon.compute_cell(c, j, i);
				c.neighbors(neigh);  // Computes a vector list of neighbors of the current particle under consideration
				c.face_vertices(vert);
				fng = 0;
				// loop over the neighbors
				for (k = 0; k < neigh.size(); k++) {

					// Skip if the neighbor information is smaller than this particle's ID, to avoid double counting. This also
					// removes faces that touch the walls, since the neighbor information is set to negative numbers for these cases.
					if (neigh[k] > id) {

						find_pos(ijk, q, neigh[k], &rcon);  // find position of the i-th neighbor

						face_dist(fng, vert, x, y, z, xn, yn, zn, c);		// {1PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
						xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;	// {1PER} urcim skutecne souradnice tohoto souseda za stenou

						rcon.compute_cell(d, ijk, q);	// compute the cell of particle
						enr = V2(c, d, x, y, z, xn, yn, zn);
						t = t + enr;
						if (enr > 0) { cit++; }		// increase number of pairs creating the energy (only if energy is > 0)

					}
					fng = fng + vert[fng] + 1;
				}
			}

		}
	}

	//return t;
	return t / cit;
}


double zet_fr(voro::container_poly &con, double R) {

	// [in]		con		container with stored generators
	// [in]		R		range

	int n0 = 0;
	double v0;

	for (int j = 0; j < con.nxyz; j++) { // loop over boxes
		for (int i = 0; i < con.co[j]; i++) { // loop over particles in considered box

		}
	}

	return 0;

}