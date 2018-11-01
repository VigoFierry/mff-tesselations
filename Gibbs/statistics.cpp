#include "Header.h"

using namespace voro;

// fce corr_coef spocte vyberovy korelacni koeficient dvojic (V_i,V_j), kde V_i a V_j jsou objemy sousednich zrn (pres vsechny dvojice sousednich zrn)
double corr_coef(voro::container &rcon)
{
	// [in]		con		the container with stored particles.

	double X, Y, XY, X2, Y2; // promenne pro ukladani prislusnych sum

	double dec = 10000000000000; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI

	unsigned int i, fng;
	int ijk, q, id;
	long n = 0;
	double a, b;
	double x, y, z, nx, ny, nz;
	voronoicell_neighbor c;  // bunka
	voronoicell d;
	std::vector<int> neigh, vert;  // vektor pro prirazeni ID sousedu, vrcholu

								   // ted chci projit pres vsechny castice, urcit jejich sousedy v mozaice, a pro ty s vyssim ID nez je samotne
								   // ID cislo aktualni castice (vyhnu se tak dvojitemu zapocitani) spocitat dvojice objemu

								   // pro prochazeni pouziji tridu loop
	c_loop_all clo(rcon);          // loop pres vsechny castice
								   // c_loop_order clo(con, po);    // loop jen pres vybrane

	X = 0; Y = 0; XY = 0; X2 = 0; Y2 = 0;  // X is sum of x_i, Y is sum of y_i, XY is sum of x_i*y_i, X2 is sum of (x_i)^2

	if (clo.start()) do if (rcon.compute_cell(c, clo)) {
		// start() - sets the class to consider the first particle 
		// compute_cell(c, clo) - return true if the computation of the Voronoi cell for a particle currently being referenced by a loop class was succesful 

		id = clo.pid();     // Get the ID of the current particle under consideration
		clo.pos(x, y, z);   // {PER} Get the position of the current particle under consideration

		c.neighbors(neigh);		// Computes a vector list of neighbors of the current particle under consideration
		c.face_vertices(vert);  // {PER} Computes list of vertices for each face of the cell

		a = c.volume();
		fng = 0;			 // {PER} 

		for (i = 0; i < neigh.size(); i++) {  // loop over the neighbors

											  // Skip if the neighbor information is smaller than this particle's ID, to avoid double counting. 
											  // # In case of non-periodicity this also removes faces that touch the walls, since the neighbor information is set to negative numbers for these cases.
											  // # In case of periodicity we need to decide which couple of neighbors to use according to its barycentrum
			if (neigh[i] > id) {

				face_dist(fng, vert, x, y, z, nx, ny, nz, c);		 // {PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
				nx = x + 2 * nx; ny = y + 2 * ny; nz = z + 2 * nz; // {PER} urcim skutecne souradnice tohoto souseda za stenou

				find_pos(ijk, q, neigh[i], &rcon);			 // najdu pozici aktualniho souseda 
				rcon.compute_cell(d, ijk, q);				 // spocte sousedni bunku
				b = d.volume();								 // a jeji objem

															 // "pres vsechny dvojice sousednich zrn" !!??? - tzn. vcetne tech periodickych nebo ne?, 
															 //      popr ridit se pravidlem teziste???? ?????????????????????????????????????????????????????????

															 // jeslize je vzd souseda od centra dvojnasobna vzd od prislusne steny (tj. nedoslo k preklopeni) 
															 // if (point_dist(x, y, z, nx, ny, nz) == 2 * dist) {

															 //round ve fci face_dist:	nx = round(nx * dec) / dec; ny = round(ny * dec) / dec; nz = round(nz * dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
															 //rounding err	
//				if (nx == rcon.p[ijk][3 * q] && ny == rcon.p[ijk][3 * q + 1] && nz == rcon.p[ijk][3 * q + 2]) {
				if (nx > 0 && nx < 1 && ny > 0 && ny < 1 && nz > 0 && nz < 1) {
					// {PER} shoduje-li se pozice souseda se skutecnou pozici za stenou (tj. nedoslo k preklopeni)
					//				if (nx > 0 && nx < 1 && ny > 0 && ny < 1 && nz > 0 && nz < 1 ) {
					// {PER} je-li pozice souseda v okne (tj. nedoslo k preklopeni)
					X += a;
					Y += b;
					XY += a*b;
					X2 += a*a;
					Y2 += b*b;

					n++;
				}
				else {
					// {PER} skutecne souradnice souseda jsou nx,ny,nz
					if (barycentrum(x, y, z, nx, ny, nz, a, b)) { // {PER} je-li teziste uvnitr okna pak maji tyto bunky vliv

						X += a;
						Y += b;
						XY += a*b;
						X2 += a*a;
						Y2 += b*b;

						n++;
					}
					//else { std::cout << nx << " " << ny << " " << nz << " \n"; } ///////////////////////////////////////
				}
			}
			fng = fng + vert[fng] + 1;  // {PER} set actual position in vector of face vertices
		}

	} while (clo.inc());    // inc() - finds the next particle to test.

							//long n = n1 + n2; ///////////////////////////////////////////////////////////////////////////////////////////////////

							//std::cout << X << " " << Y << " " << XY << " " << X2 << " " << Y2 << " \n"; /////////////////////////////////////////
							//std::cout << n1 << " " << n2 << " \n";  /////////////////////////////////////////////////////////////////////////////
							//std::cout << XY - X*Y / n << " " << X2 - X*X / n << " " << Y2 - Y*Y / n << " \n";   /////////////////////////////////

	x = round((XY - X*Y / n) *dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	y = round((X2 - X*X / n)*dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	z = round((Y2 - Y*Y / n)*dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI

	if (y == 0) { std::cout << "Sxy = " << x << "; Sz = " << z << " ; Sy = 0 "; return 0; }
	if (z == 0) { std::cout << "Sxy = " << x << "; Sy = " << y << " ; Sz = 0 "; return 0; }
	return x / sqrt(y*z);
}

double corr_coef(voro::container_poly &rcon)
{
	// [in]		con		the container with stored particles.

	double X, Y, XY, X2, Y2; // promenne pro ukladani prislusnych sum

	double dec = 10000000000000; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI

	unsigned int i, fng;
	int ijk, q, id;
	long n = 0;
	double a, b;
	double x, y, z, nx, ny, nz;
	voronoicell_neighbor c;  // bunka
	voronoicell d;
	std::vector<int> neigh, vert;  // vektor pro prirazeni ID sousedu, vrcholu

								   // ted chci projit pres vsechny castice, urcit jejich sousedy v mozaice, a pro ty s vyssim ID nez je samotne
								   // ID cislo aktualni castice (vyhnu se tak dvojitemu zapocitani) spocitat dvojice objemu

								   // pro prochazeni pouziji tridu loop
	c_loop_all clo(rcon);          // loop pres vsechny castice
								   // c_loop_order clo(con, po);    // loop jen pres vybrane

	X = 0; Y = 0; XY = 0; X2 = 0; Y2 = 0;  // X is sum of x_i, Y is sum of y_i, XY is sum of x_i*y_i, X2 is sum of (x_i)^2

	if (clo.start()) do if (rcon.compute_cell(c, clo)) {
		// start() - sets the class to consider the first particle 
		// compute_cell(c, clo) - return true if the computation of the Voronoi cell for a particle currently being referenced by a loop class was succesful 

		id = clo.pid();     // Get the ID of the current particle under consideration
		clo.pos(x, y, z);   // {PER} Get the position of the current particle under consideration

		c.neighbors(neigh);		// Computes a vector list of neighbors of the current particle under consideration
		c.face_vertices(vert);  // {PER} Computes list of vertices for each face of the cell

		a = c.volume();
		fng = 0;			 // {PER} 

		for (i = 0; i < neigh.size(); i++) {  // loop over the neighbors

											  // Skip if the neighbor information is smaller than this particle's ID, to avoid double counting. 
											  // # In case of non-periodicity this also removes faces that touch the walls, since the neighbor information is set to negative numbers for these cases.
											  // # In case of periodicity we need to decide which couple of neighbors to use according to its barycentrum
			if (neigh[i] > id) {

				face_dist(fng, vert, x, y, z, nx, ny, nz, c);		 // {PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
				nx = x + 2 * nx; ny = y + 2 * ny; nz = z + 2 * nz; // {PER} urcim skutecne souradnice tohoto souseda za stenou

				find_pos(ijk, q, neigh[i], &rcon);			 // najdu pozici aktualniho souseda 
				rcon.compute_cell(d, ijk, q);				 // spocte sousedni bunku
				b = d.volume();								 // a jeji objem

															 // "pres vsechny dvojice sousednich zrn" !!??? - tzn. vcetne tech periodickych nebo ne?, 
															 //      popr ridit se pravidlem teziste???? ?????????????????????????????????????????????????????????

															 // jeslize je vzd souseda od centra dvojnasobna vzd od prislusne steny (tj. nedoslo k preklopeni) 
															 // if (point_dist(x, y, z, nx, ny, nz) == 2 * dist) {

															 //round ve fci face_dist:	nx = round(nx * dec) / dec; ny = round(ny * dec) / dec; nz = round(nz * dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
															 //rounding err	
															 //				if (nx == rcon.p[ijk][3 * q] && ny == rcon.p[ijk][3 * q + 1] && nz == rcon.p[ijk][3 * q + 2]) {
				if (nx > 0 && nx < 1 && ny > 0 && ny < 1 && nz > 0 && nz < 1) {
					// {PER} shoduje-li se pozice souseda se skutecnou pozici za stenou (tj. nedoslo k preklopeni)
					//				if (nx > 0 && nx < 1 && ny > 0 && ny < 1 && nz > 0 && nz < 1 ) {
					// {PER} je-li pozice souseda v okne (tj. nedoslo k preklopeni)
					X += a;
					Y += b;
					XY += a*b;
					X2 += a*a;
					Y2 += b*b;

					n++;
				}
				else {
					// {PER} skutecne souradnice souseda jsou nx,ny,nz
					if (barycentrum(x, y, z, nx, ny, nz, a, b)) { // {PER} je-li teziste uvnitr okna pak maji tyto bunky vliv

						X += a;
						Y += b;
						XY += a*b;
						X2 += a*a;
						Y2 += b*b;

						n++;
					}
					//else { std::cout << nx << " " << ny << " " << nz << " \n"; } ///////////////////////////////////////
				}
			}
			fng = fng + vert[fng] + 1;  // {PER} set actual position in vector of face vertices
		}

	} while (clo.inc());    // inc() - finds the next particle to test.

							//long n = n1 + n2; ///////////////////////////////////////////////////////////////////////////////////////////////////

							//std::cout << X << " " << Y << " " << XY << " " << X2 << " " << Y2 << " \n"; /////////////////////////////////////////
							//std::cout << n1 << " " << n2 << " \n";  /////////////////////////////////////////////////////////////////////////////
							//std::cout << XY - X*Y / n << " " << X2 - X*X / n << " " << Y2 - Y*Y / n << " \n";   /////////////////////////////////

	x = round((XY - X*Y / n) *dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	y = round((X2 - X*X / n)*dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	z = round((Y2 - Y*Y / n)*dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI

	if (y == 0) { std::cout << "Sxy = " << x << "; Sz = " << z << " ; Sy = 0 "; return 0; }
	if (z == 0) { std::cout << "Sxy = " << x << "; Sy = " << y << " ; Sz = 0 "; return 0; }
	return x / sqrt(y*z);
}


void con_stats(voro::container &rcon)
{
	// [in]		con				the container with stored particles

	FILE *f;
	f = fopen("../data/con_stats.txt", "w");

	if (f == NULL) { std::cout << "NELZE zapsat con_stats \n"; }

	rcon.total_particles();		// pocet bunek
	// pocet sten
	// pocet hran
	// pocet vrcholu
}

void con_stats(voro::container_poly &con) {
	int i, j, n;
	int citac = 0;
	double rad = 0;

	//FILE *f;
	//f = fopen("datacon.txt", "w");
	//f = fopen("../data/d200_1.txt", "w");
	//if (f == NULL) { std::cout << "NELZE zapsat container \n"; }

	for (j = 0; j < con.nxyz; j++) {
		for (i = 0; i < con.co[j]; i++) {
			if (con.p[j][4 * i + 3] == 0) { citac++; }
			rad = rad + con.p[j][4 * i + 3];
		}
	}
	n = con.total_particles();

	std::cout << "Pocet nulovych polomeru: " << citac << "\n";
	std::cout << "Prumerny polomer: " << rad / n << "\n";
}

void cell_stats(voro::container &rcon)
{
	// [in]		con				the container with stored particles

	int i, j;
	voronoicell_neighbor c;
	// soubor pro zapisovani dat:
	FILE *f;
	//C:\Users\Vigo\Documents\Visual Studio 2015\Projects\Gibbs - Voronoi\Gibbs - Voronoi
	//f = fopen("C://Users/Vigo/Dokumenty/Visual Studio 2015/Projects/Gibbs-Voronoi/Gibbs-Voronoi/cell_stats.txt", "w");
	f = fopen("cell_stats.txt", "w");
		
	//f = fopen("../data/c200_1_1.txt", "w");
	if (f == NULL) { std::cout << "NELZE zapsat cell_stats \n"; }
   
	// loop pres vsechny castice:
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes

		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			rcon.compute_cell(c, j, i);		// spocte se bunka aktualne uvazovane castice

											// fprintf(f, "%d %g %g %g \n", rcon.id[j][i], con.p[j][3 * i], con.p[j][3 * i + 1], con.p[j][3 * i + 2]);

			fprintf(f, "%d %g %g %g %d %d %g \n",
				rcon.id[j][i],				// ID
				// bunka:
				c.volume(),					// objem
				c.surface_area(),			// povrch
				c.total_edge_distance(),	// delka hran
				c.number_of_faces(),		// pocet sten = pocet sousedu
				c.number_of_edges(),		// pocet hran
				// pocet vrcholu bunky = pocet hran - pocet sten + 2
				c.max_radius_squared()		// druha mocnina maximalni vzdalenosti vrcholu od stredu
			);
		}
	}
	fclose(f);
}

void cell_stats(voro::container_poly &rcon)
{
	// [in]		con				the container with stored particles

	int i, j;
	bool grain;
	voronoicell_neighbor c;
	// soubor pro zapisovani dat:
	FILE *f;
	//C:\Users\Vigo\Documents\Visual Studio 2015\Projects\Gibbs - Voronoi\Gibbs - Voronoi
	//f = fopen("C://Users/Vigo/Dokumenty/Visual Studio 2015/Projects/Gibbs-Voronoi/Gibbs-Voronoi/cell_stats.txt", "w");
	f = fopen("cell_statsHIST.txt", "w");

	//f = fopen("../data/c200_1_1.txt", "w");
	if (f == NULL) { std::cout << "NELZE zapsat cell_stats \n"; }

	// loop pres vsechny castice:
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes

		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			grain = rcon.compute_cell(c, j, i);		// spocte se bunka aktualne uvazovane castice

											// fprintf(f, "%d %g %g %g \n", rcon.id[j][i], con.p[j][3 * i], con.p[j][3 * i + 1], con.p[j][3 * i + 2]);
			if (grain == true) {
				fprintf(f, "%d %g %g %g %d %d %g \n",
					rcon.id[j][i],				// ID
												// bunka:
					c.volume(),					// objem
					c.surface_area(),			// povrch
					c.total_edge_distance(),	// delka hran
					c.number_of_faces(),		// pocet sten = pocet sousedu
					c.number_of_edges(),		// pocet hran
												// pocet vrcholu bunky = pocet hran - pocet sten + 2
					c.max_radius_squared()		// druha mocnina maximalni vzdalenosti vrcholu od stredu
				);
			}
		}
	}
	fclose(f);
}

void neigh_list(voro::container_poly &rcon)
{
	// [in]		con				the container with stored particles

	int i, j,k;
	bool grain;
	voronoicell_neighbor c;

	std::vector<int> neigh;
	// soubor pro zapisovani dat:
	FILE *f;
	//C:\Users\Vigo\Documents\Visual Studio 2015\Projects\Gibbs - Voronoi\Gibbs - Voronoi
	//f = fopen("C://Users/Vigo/Dokumenty/Visual Studio 2015/Projects/Gibbs-Voronoi/Gibbs-Voronoi/cell_stats.txt", "w");
	f = fopen("Nlist.txt", "w");

	//f = fopen("../data/c200_1_1.txt", "w");
	if (f == NULL) { std::cout << "NELZE zapsat Nlist \n"; }

	// loop pres vsechny castice:
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes

		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			grain = rcon.compute_cell(c, j, i);		// spocte se bunka aktualne uvazovane castice
			if (grain == true) {
				c.neighbors(neigh);
				fprintf(f, "%d ", rcon.id[j][i]);
				fprintf(f, "%d ", neigh.size());
				// fprintf(f, "%d %g %g %g \n", rcon.id[j][i], con.p[j][3 * i], con.p[j][3 * i + 1], con.p[j][3 * i + 2]);
				for (k = 0; k < neigh.size(); k++) {

					fprintf(f, "%d ",
						neigh[k]
					);
				}
				fprintf(f, " \n");
			}
		}
	}
	fclose(f);
}

void face_stats(voro::container &rcon)
{
	// [in]		con				the container with stored particles

	int i, j, ijk, q;
	unsigned int k;
	double vol;
	voronoicell_neighbor c;
	voronoicell d;
	std::vector<int> neigh, ford;
	std::vector<double> fare, fper;
	// soubor pro zapisovani dat:
	FILE *f;
	f = fopen("face_stats.txt", "w");
		
	//f = fopen("../data/f200_1_1.txt", "w");
	if (f == NULL) { std::cout << "NELZE zapsat face_stats \n"; }
	

	// loop pres vsechny castice:
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes

		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			rcon.compute_cell(c, j, i);		// spocte se bunka aktualne uvazovane castice
			c.neighbors(neigh);
			c.face_areas(fare); c.face_orders(ford); c.face_perimeters(fper);
			vol = c.volume();

			for (k = 0; k < neigh.size(); k++) {	// loop over neighbors = loop over faces
				if (rcon.id[j][i] < neigh[k]) {		// if ID of neighbor is higher than id of current cell continue (to avoid duble counting)
					find_pos(ijk, q, neigh[k], &rcon);
					rcon.compute_cell(d, ijk, q);
					fprintf(f, "%g %d %g %g %g\n",
						fare[k],		// povrch steny
						ford[k],		// pocet vrcholu steny = pocet hran steny
						fper[k],		// obvod steny
						vol,			// objem prvni bunky
						d.volume()		// objem druhe bunky
					);
				}
			}
		}
	}
	fclose(f);
}

void face_stats(voro::container_poly &rcon)
{
	// [in]		con				the container with stored particles

	int i, j, ijk, q;
	bool grain, frain;
	unsigned int k;
	double vol;
	voronoicell_neighbor c;
	voronoicell d;
	std::vector<int> neigh, ford;
	std::vector<double> fare, fper;
	// soubor pro zapisovani dat:
	FILE *f;
	f = fopen("face_stats.txt", "w");

	//f = fopen("../data/f200_1_1.txt", "w");
	if (f == NULL) { std::cout << "NELZE zapsat face_stats \n"; }


	// loop pres vsechny castice:
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes

		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			grain = rcon.compute_cell(c, j, i);		// spocte se bunka aktualne uvazovane castice
			if (grain  == true) {
				c.neighbors(neigh);
				c.face_areas(fare); c.face_orders(ford); c.face_perimeters(fper);
				vol = c.volume();

				for (k = 0; k < neigh.size(); k++) {	// loop over neighbors = loop over faces
					if (rcon.id[j][i] < neigh[k]) {		// if ID of neighbor is higher than id of current cell continue (to avoid duble counting)
						find_pos(ijk, q, neigh[k], &rcon);
						frain = rcon.compute_cell(d, ijk, q);
						if (frain = true) {
							fprintf(f, "%g %d %g %g %g\n",
								fare[k],		// povrch steny
								ford[k],		// pocet vrcholu steny = pocet hran steny
								fper[k],		// obvod steny
								vol,			// objem prvni bunky
								d.volume()		// objem druhe bunky
							);
						}
					}
				}
			}
		}
	}
	fclose(f);
}

void edge_stats(voro::container &rcon)
{
	// [in]		con				the container with stored particles

	// edge is determined by three neighboring cells
	int i, j, k, l;
	int ijk, q, p, i1, i2;
	voronoicell_neighbor c,d1,d2;
	std::vector<int> neigh, neigh1, neigh2;
	std::vector<double> lengths, angles, angles1, angles2, fare, fare1;
	// soubor pro zapisovani dat:
	FILE *f;
	f = fopen("edge_stats.txt", "a");
	//f = fopen("../data/edge_stats.txt", "w");

	if (f == NULL) { std::cout << "NELZE zapsat edge_stats \n"; }

	// loop pres vsechny castice:
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes

		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			rcon.compute_cell(c, j, i);		// spocte se bunka aktualne uvazovane castice
			c.neighbors(neigh);
			c.face_areas(fare);
			edge_lengths(lengths, c);
			dihedral_angles(angles, c);
			// to avoid triple counting we need IF conditions
			for (k = 0; k < neigh.size(); k++) {	// loop over neighbors = loop over faces
				if (rcon.id[j][i] < neigh[k]) {				// the first IF 
					find_pos(ijk, q, neigh[k], &rcon);
					rcon.compute_cell(d1, ijk, q);
					d1.neighbors(neigh1);
					d1.face_areas(fare1);
					dihedral_angles(angles1, d1);
					for (l = 0; l < neigh.size(); l++) {	// loop over neighbors again = loop over faces again
						if (neigh[k] < neigh[l]) {					// the second IF 
							find_pos(ijk, q, neigh[l], &rcon);
							p = common_edge(c, k, l);
							if (p > -1) {					// maji-li spolecnou hranu (trojice sousedicich bunek nemusi urcovat hranu) 
								// predpokladame-li ze hrana je urcena 3 bunkami/stenami, tak je toto zbytecne (podminka s common edge)
								if (are_neighbors(d1, ijk, q, &rcon)) {		// najdeme tak trojici sousedicich bunek urcujicich hranu (techto bunek ale muze byt vice nez 3 !!!)
											// navic diky IF podminkam se kazda hrana vyskytne jen jednou
									// prvni bunku pouziju, u druhe a treti mi staci vedet, ktere steny jim prislusi
									// tyto dve steny urcuji take hranu, staci urcit spolecne vrcholy a poradi teto hrany ve strukture ed -> p
									// p = common_edge(c, k, l);
									fprintf(f, "%g ", lengths[p]); // delka hrany

									// PRUSER

									// dale me zajimaji povrchy vsech tri sten vybihajicich z teto hrany, nato potrebuji dve bunky
									fprintf(f, "%g %g ", fare[k], fare[l]); // obsahy sten, ktere lze najit v prvni bunce
									// s dihedralnimi uhly je to slozitejsi, prtz kazdy z nich prislusi jine bunce
									// kdybych nepouzil podminky IF spocetl bych vsechny dihedralni uhly snaze, ale nedokazal bych rozlisit, kterym
									// hranam patri
									// vytisknu az pozdeji - fprintf(f, "%g ", angles[p]);  // dihedralni uhel prislusejici prvni bunce c

									// pro druhou bunku: najdu zbyle dve bunky mezi jejimi sousedy, tim urcim steny, pak najdu hranu atd
									i1 = find_index(neigh1, rcon.id[j][i]);
									i2 = find_index(neigh1, neigh[l]);
									fprintf(f, "%g ", fare1[i2]); // obsah treti steny
									fprintf(f, "%g ", angles[p]);  // dihedralni uhel prislusejici prvni bunce c
									p = common_edge(d1, i1, i2);
									fprintf(f, "%g ", angles1[p]);  // dihedralni uhel prislusejici druhe bunce d1

									// pro treti bunku analogicky:
									rcon.compute_cell(d2, ijk, q);
									d2.neighbors(neigh2);
									dihedral_angles(angles2, d2);
									i1 = find_index(neigh2, rcon.id[j][i]);
									i2 = find_index(neigh2, neigh[k]);
									p = common_edge(d2, i1, i2);
									fprintf(f, "%g \n", angles2[p]);  // dihedralni uhel prislusejici treti bunce d2
								}
								else {		// nevime kolik bunek obsahuje tuto hranu, vypiseme jen udaje z prvni bunky
									// fprintf(f, "%g ", lengths[p]);			// delka hrany
									// fprintf(f, "%g %g ", fare[k], fare[l]); // obsahy sten, ktere lze najit v prvni bunce
									// fprintf(f, "NA ");			//  nezname kolik dalsich sten obsahuje tutez hranu
									// fprintf(f, "%g ", angles[p]);			 // dihedralni uhel prislusejici prvni bunce c
									// fprintf(f, "NA NA NA \n");
									fprintf(f, "NA \n");
								}
							}
						}//END (if; k,l)
					}//END (for; l)
				}//END (if; i,j,k)
			}//END (for; k)
		}//END (for; i)
	}//END (for; j)
	fclose(f);
}

void edge_stats(voro::container_poly &rcon)
{
	// [in]		con				the container with stored particles

	// edge is determined by three neighboring cells
	int i, j, k, l;
	int ijk, q, p, i1, i2;
	bool grain;
	voronoicell_neighbor c, d1, d2;
	std::vector<int> neigh, neigh1, neigh2;
	std::vector<double> lengths, angles, angles1, angles2, fare, fare1;
	// soubor pro zapisovani dat:
	FILE *f;
	f = fopen("edge_stats.txt", "a");
	//f = fopen("../data/edge_stats.txt", "w");

	if (f == NULL) { std::cout << "NELZE zapsat edge_stats \n"; }

	// loop pres vsechny castice:
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes

		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			grain = rcon.compute_cell(c, j, i);		// spocte se bunka aktualne uvazovane castice
			if (grain == true) {
				c.neighbors(neigh);
				c.face_areas(fare);
				edge_lengths(lengths, c);
				dihedral_angles(angles, c);
				// to avoid triple counting we need IF conditions
				for (k = 0; k < neigh.size(); k++) {	// loop over neighbors = loop over faces
					if (rcon.id[j][i] < neigh[k]) {				// the first IF 
						find_pos(ijk, q, neigh[k], &rcon);
						rcon.compute_cell(d1, ijk, q);
						d1.neighbors(neigh1);
						d1.face_areas(fare1);
						dihedral_angles(angles1, d1);
						for (l = 0; l < neigh.size(); l++) {	// loop over neighbors again = loop over faces again
							if (neigh[k] < neigh[l]) {					// the second IF 
								find_pos(ijk, q, neigh[l], &rcon);
								p = common_edge(c, k, l);
								if (p > -1) {					// maji-li spolecnou hranu (trojice sousedicich bunek nemusi urcovat hranu) 
																// podminka p>-1, tj podminka s common_edge je dulezitejsi nez podminka are_neighbors (ta je redundantni)
									if (are_neighbors(d1, ijk, q, &rcon)) {		// najdeme tak trojici sousedicich bunek urcujicich hranu (techto bunek ale muze byt vice nez 3 !!!)
																				// navic diky IF podminkam se kazda hrana vyskytne jen jednou
																				// prvni bunku pouziju, u druhe a treti mi staci vedet, ktere steny jim prislusi
																				// tyto dve steny urcuji take hranu, staci urcit spolecne vrcholy a poradi teto hrany ve strukture ed -> p
																				// p = common_edge(c, k, l);
										fprintf(f, "%g ", lengths[p]); // delka hrany

																	   // PRUSER

																	   // dale me zajimaji povrchy vsech tri sten vybihajicich z teto hrany, nato potrebuji dve bunky
										fprintf(f, "%g %g ", fare[k], fare[l]); // obsahy sten, ktere lze najit v prvni bunce
																				// s dihedralnimi uhly je to slozitejsi, prtz kazdy z nich prislusi jine bunce
																				// kdybych nepouzil podminky IF spocetl bych vsechny dihedralni uhly snaze, ale nedokazal bych rozlisit, kterym
																				// hranam patri
																				// vytisknu az pozdeji - fprintf(f, "%g ", angles[p]);  // dihedralni uhel prislusejici prvni bunce c

																				// pro druhou bunku: najdu zbyle dve bunky mezi jejimi sousedy, tim urcim steny, pak najdu hranu atd
										i1 = find_index(neigh1, rcon.id[j][i]);
										if (i1 == -1) { std::cout << "ERROR: edge statistics - find index \n"; return; }
										i2 = find_index(neigh1, neigh[l]);
										if (i2 == -1) { std::cout << "ERROR: edge statistics - find index \n"; return; }
										fprintf(f, "%g ", fare1[i2]); // obsah treti steny
										fprintf(f, "%g ", angles[p]);  // dihedralni uhel prislusejici prvni bunce c
										p = common_edge(d1, i1, i2);
										if (p == -1) { std::cout << "ERROR: edge statistics - common edge \n"; return; }
										fprintf(f, "%g ", angles1[p]);  // dihedralni uhel prislusejici druhe bunce d1

																		// pro treti bunku analogicky:
										rcon.compute_cell(d2, ijk, q);
										d2.neighbors(neigh2);
										dihedral_angles(angles2, d2);
										i1 = find_index(neigh2, rcon.id[j][i]);
										if (i1 == -1) { std::cout << "ERROR: edge statistics - find index \n"; return; }
										i2 = find_index(neigh2, neigh[k]);
										if (i1 == -1) { std::cout << "ERROR: edge statistics - find index \n"; return; }
										p = common_edge(d2, i1, i2);
										if (p == -1) { std::cout << "ERROR: edge statistics - common edge \n"; return; }
										fprintf(f, "%g \n", angles2[p]);  // dihedralni uhel prislusejici treti bunce d2
									}
									else {		// nevime kolik bunek obsahuje tuto hranu, vypiseme jen udaje z prvni bunky
												// fprintf(f, "%g ", lengths[p]);			// delka hrany
												// fprintf(f, "%g %g ", fare[k], fare[l]); // obsahy sten, ktere lze najit v prvni bunce
												// fprintf(f, "NA ");			//  nezname kolik dalsich sten obsahuje tutez hranu
												// fprintf(f, "%g ", angles[p]);			 // dihedralni uhel prislusejici prvni bunce c
												// fprintf(f, "NA NA NA \n");
										fprintf(f, "NA \n");
									}
								}
							}//END (if; k,l)
						}//END (for; l)
					}//END (if; i,j,k)
				}//END (for; k)
			} //END (if; grain)
		}//END (for; i)
	}//END (for; j)
	fclose(f);
}

void vertex_stats(voro::container &rcon)	// NEODLADENA
{
	// [in]		con				the container with stored particles

	// vertex is determined by four cells
	// vrchol je ale soucasti tri sten a kazda z nich urcuje jednu sousedni bunku --> tak dostanu 4 bunky
	int i, j, k, l, m;
	int ijk, q;
	int nv, nf, n, p, f1, f2, f3, n1, n2, n3, p1, p2, p3, v1, v2;
	voronoicell_neighbor c, d;
	std::vector<int> neigh, fvert, ford, vord, faces, neigh1;
	std::vector<double> fle, len;
	// soubor pro zapisovani dat:
	FILE *f;
	//f = fopen("vertex_stats.txt", "w");
	f = fopen("../data/vertex_stats.txt", "w");

	if (f == NULL) { std::cout << "NELZE zapsat vertex_stats \n"; }

	// loop pres vsechny castice:
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes

		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			rcon.compute_cell(c, j, i);		// spocte se bunka aktualne uvazovane castice
			c.neighbors(neigh);
			c.face_vertices(fvert);
			c.face_orders(ford);
			c.vertex_orders(vord);
			face_edge_lengths(fle, c);
			nf = c.number_of_faces();
			nv = c.number_of_edges() - nf + 2; // nv ... number of vertices

			for (k = 0; k < nv; k++) {		// loop over vertices
				p = 1; n = 0;
				for (l = 0; l < nf; l++) {		// find three faces of the actual vertex
					for (m = 0; m < ford[l]; m++) {
						if(k == fvert[p+m]){
							// l1,l2,l3 ... delka hrany vybihajici z vrcholu (!!! bude fungovat jen tehdy, pokud k prochazeni
							// hran v kazde stene dochazi ve stejnem smeru (po nebo proti smeru hodinovych rucicek) !!!!!!!!!!!
							/*
							if (n == 1) { f1 = l; n++; l1 = fle[p + m]; }		// the first face
							if (n == 2) { f2 = l; n++; l2 = fle[p + m]; }		// the second face
							if (n == 3) { f3 = l; n++; l3 = fle[p + m]; }		// the third face
							if (n == 4){ std::cout << "VERTEX: the FOURTH FACE FOUND! \n"; std::cin.ignore(); }	// PRUSER
							*/
							faces[n] = l;
							len[n] = fle[p + m];

						}
					}
					p = p + ford[l];
				}
				n1 = neigh[faces[0]];
				for (l = 0; l < vord[k]; l++) {		// n1 nejmensi index souseda
					if (neigh[faces[l]] < n1) { n1 = neigh[faces[l]]; }
				}
				// n1 = neigh[f1]; n2 = neigh[f2]; n3 = neigh[f3];		// determine neighbors
				// min_max(n1, n2, n3);		// n1 < n2 < n3

				if (rcon.id[j][i] < n1) {	// to avoid quatro counting continue only with the cell with the smallest ID
					// fprintf(f, "%g %g %g ", l1, l2, l3);
					fprintf(f, "%d ", vord[k]);			// rad vrcholu
					for (l = 0; l < vord[k]; l++) {
						fprintf(f, "%g ", len[l]);		// delky hran patricich do prvni bunky
					}
					// zbyva urcit delku ctvrte hrany, ktera se nachazi v jine bunce
					find_pos(ijk, q, n1, &rcon);
					rcon.compute_cell(d, ijk, q);
					d.neighbors(neigh1);
					// 4) najdu ty tri steny v teto bunce
					f1 = find_index(neigh1, rcon.id[j][i]);
					if (f1 == -1) { std::cout << "ERROR: Stena nenalezena (f1)! \n"; break; }
					f2 = find_index(neigh1, n2);	// tyto steny budou mit opet jeden spolecny vrchol
					if (f2 == -1) { std::cout << "ERROR: Stena nenalezena (f2)! \n"; break; }
					f3 = find_index(neigh1, n3);
					if (f3 == -1) { std::cout << "ERROR: Stena nenalezena (f3)! \n"; break; }
					// 5) najdu spolecny vrchol techto sten
					d.face_orders(ford);
					p1 = 0;
					for (l = 0; l < f1; l++) {
						p1 = p1 + ford[l];
					}
					p1 = p1 + f1 + 1;			// p1 - prvni vrchol steny f1
					p2 = 0;
					for (l = 0; l < f2; l++) {
						p2 = p2 + ford[l];
					}
					p2 = p2 + f2 + 1;			// p2 - prvni vrchol steny f2
					p3 = 0;
					for (l = 0; l < f3; l++) {
						p3 = p3 + ford[l];
					}
					p3 = p3 + f3 + 1;			// p3 - prvni vrchol steny f3

					d.face_vertices(fvert);
					for (n = 0; n < ford[f1]; n++) {				//najdi spolecny vrchol techto sten
						for (l = 0; l < ford[f2]; l++) {
							for (m = 0; m < ford[f3]; m++) {
								if (fvert[p1 + k] == fvert[p2 + l] && fvert[p2 + l] == fvert[p3 + m]) { v1 = fvert[p1 + k]; break; }
							}
						}
					}
					// 6) Ze tri sousedu tohoto vrcholu budou dva lezet v prvni stene a treti nikoliv. Ten najdu.
					for (l = 0; l < d.nu[v1]; l++) {			// najdi sousedni vrcchol, ktery nelezi v prvni stene
						v2 = d.ed[v1][l];
						n = 0;
						for (m = 0; m < ford[f1]; m++) {
							if (v2 == fvert[p1 + l]) { n = 1; }  // lezi v prvni stene, pokracuj
						}
						if (n == 0) { break; }					// nasli jsme vrchol, ktery neni v prvni stene, ukonci cyklus
					}

					// 7) Pro dva nalezene vrcholy, tj. pro posledni hranu, urci jeji poradi ve strukture ed a spocti jeji delku.
					min_max(v1, v2);			// v1 >= v2

					m = 0;										// najdi poradi teto hrany
					for (l = 0; l < v2; l++) {
						m = m + d.nu[l];
					}
					l = 0;
					m++;
					while (d.ed[v2][l] != v1) {
						m++; l++; if (l > d.nu[v2]-1) { break; std::cout << " ERROR: Edge not found! \n"; }
					}
					edge_lengths(fle, d);
					fprintf(f, "%g \n", fle[m]);	// delka posledni hrany

				}
			}

		}
	}
	fclose(f);
}

void vertex_stats(voro::container_poly &rcon)	// NEODLADENA
{
	// [in]		con				the container with stored particles

	// vertex is determined by four cells
	// vrchol je ale soucasti tri sten a kazda z nich urcuje jednu sousedni bunku --> tak dostanu 4 bunky
	int i, j, k, l, m;
	int ijk, q;
	int nv, nf, n, p, f1, f2, f3, n1, n2, n3, p1, p2, p3, v1, v2;
	bool grain;
	voronoicell_neighbor c, d;
	std::vector<int> neigh, fvert, ford, vord, faces, neigh1;
	std::vector<double> fle, len;
	// soubor pro zapisovani dat:
	FILE *f;
	//f = fopen("vertex_stats.txt", "w");
	f = fopen("../data/vertex_stats.txt", "w");

	if (f == NULL) { std::cout << "NELZE zapsat vertex_stats \n"; }

	// loop pres vsechny castice:
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes

		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			rcon.compute_cell(c, j, i);		// spocte se bunka aktualne uvazovane castice
			c.neighbors(neigh);
			c.face_vertices(fvert);
			c.face_orders(ford);
			c.vertex_orders(vord);
			face_edge_lengths(fle, c);
			nf = c.number_of_faces();
			nv = c.number_of_edges() - nf + 2; // nv ... number of vertices

			for (k = 0; k < nv; k++) {		// loop over vertices
				p = 1; n = 0;
				for (l = 0; l < nf; l++) {		// find three faces of the actual vertex
					for (m = 0; m < ford[l]; m++) {
						if (k == fvert[p + m]) {
							// l1,l2,l3 ... delka hrany vybihajici z vrcholu (!!! bude fungovat jen tehdy, pokud k prochazeni
							// hran v kazde stene dochazi ve stejnem smeru (po nebo proti smeru hodinovych rucicek) !!!!!!!!!!!
							/*
							if (n == 1) { f1 = l; n++; l1 = fle[p + m]; }		// the first face
							if (n == 2) { f2 = l; n++; l2 = fle[p + m]; }		// the second face
							if (n == 3) { f3 = l; n++; l3 = fle[p + m]; }		// the third face
							if (n == 4){ std::cout << "VERTEX: the FOURTH FACE FOUND! \n"; std::cin.ignore(); }	// PRUSER
							*/
							faces[n] = l;
							len[n] = fle[p + m];

						}
					}
					p = p + ford[l];
				}
				n1 = neigh[faces[0]];
				for (l = 0; l < vord[k]; l++) {		// n1 nejmensi index souseda
					if (neigh[faces[l]] < n1) { n1 = neigh[faces[l]]; }
				}
				// n1 = neigh[f1]; n2 = neigh[f2]; n3 = neigh[f3];		// determine neighbors
				// min_max(n1, n2, n3);		// n1 < n2 < n3

				if (rcon.id[j][i] < n1) {	// to avoid quatro counting continue only with the cell with the smallest ID
											// fprintf(f, "%g %g %g ", l1, l2, l3);
					fprintf(f, "%d ", vord[k]);			// rad vrcholu
					for (l = 0; l < vord[k]; l++) {
						fprintf(f, "%g ", len[l]);		// delky hran patricich do prvni bunky
					}
					// zbyva urcit delku ctvrte hrany, ktera se nachazi v jine bunce
					find_pos(ijk, q, n1, &rcon);
					rcon.compute_cell(d, ijk, q);
					d.neighbors(neigh1);
					// 4) najdu ty tri steny v teto bunce
					f1 = find_index(neigh1, rcon.id[j][i]);
					if (f1 == -1) { std::cout << "ERROR: Stena nenalezena (f1)! \n"; break; }
					f2 = find_index(neigh1, n2);	// tyto steny budou mit opet jeden spolecny vrchol
					if (f2 == -1) { std::cout << "ERROR: Stena nenalezena (f2)! \n"; break; }
					f3 = find_index(neigh1, n3);
					if (f3 == -1) { std::cout << "ERROR: Stena nenalezena (f3)! \n"; break; }
					// 5) najdu spolecny vrchol techto sten
					d.face_orders(ford);
					p1 = 0;
					for (l = 0; l < f1; l++) {
						p1 = p1 + ford[l];
					}
					p1 = p1 + f1 + 1;			// p1 - prvni vrchol steny f1
					p2 = 0;
					for (l = 0; l < f2; l++) {
						p2 = p2 + ford[l];
					}
					p2 = p2 + f2 + 1;			// p2 - prvni vrchol steny f2
					p3 = 0;
					for (l = 0; l < f3; l++) {
						p3 = p3 + ford[l];
					}
					p3 = p3 + f3 + 1;			// p3 - prvni vrchol steny f3

					d.face_vertices(fvert);
					for (n = 0; n < ford[f1]; n++) {				//najdi spolecny vrchol techto sten
						for (l = 0; l < ford[f2]; l++) {
							for (m = 0; m < ford[f3]; m++) {
								if (fvert[p1 + k] == fvert[p2 + l] && fvert[p2 + l] == fvert[p3 + m]) { v1 = fvert[p1 + k]; break; }
							}
						}
					}
					// 6) Ze tri sousedu tohoto vrcholu budou dva lezet v prvni stene a treti nikoliv. Ten najdu.
					for (l = 0; l < d.nu[v1]; l++) {			// najdi sousedni vrcchol, ktery nelezi v prvni stene
						v2 = d.ed[v1][l];
						n = 0;
						for (m = 0; m < ford[f1]; m++) {
							if (v2 == fvert[p1 + l]) { n = 1; }  // lezi v prvni stene, pokracuj
						}
						if (n == 0) { break; }					// nasli jsme vrchol, ktery neni v prvni stene, ukonci cyklus
					}

					// 7) Pro dva nalezene vrcholy, tj. pro posledni hranu, urci jeji poradi ve strukture ed a spocti jeji delku.
					min_max(v1, v2);			// v1 >= v2

					m = 0;										// najdi poradi teto hrany
					for (l = 0; l < v2; l++) {
						m = m + d.nu[l];
					}
					l = 0;
					m++;
					while (d.ed[v2][l] != v1) {
						m++; l++; if (l > d.nu[v2] - 1) { break; std::cout << " ERROR: Edge not found! \n"; }
					}
					edge_lengths(fle, d);
					fprintf(f, "%g \n", fle[m]);	// delka posledni hrany

				}
			}

		}
	}
	fclose(f);
}

// fce stats_output by mohla byt kombinaci fci cell_stats, face_stats, edge_stats a vertex_stats a tisknout veskerou moznou informaci
// -- je to mozne vse zkombinovat do jedne fce ??? -- ano, ale musi se tisknout do nekolika souboru (zvlast soubor pro cell-,face-,edge- i 
//	vertex-stats; tj 4 soubory)

void stats_output(voro::container &rcon) // NEODLADENA
{
	// [in]		con				the container with stored particles

	int i, j, k, l;
	int ijk, q;
	int p, i1, i2;
	double v;
	voronoicell_neighbor c, d, e;
	std::vector<int> neigh, neigh1, neigh2, ford, fvert, vord;
	std::vector<double> fare, fare1, fper, vert, norm, lengths, angles, angles1, angles2;


	// soubory pro zapisovani dat:
	FILE *f; FILE *ff; FILE *fff; //FILE *ffff;
	f = fopen("cell_stats.txt", "w");
	ff = fopen("face_stats.txt", "w");
	fff = fopen("edge_stats.txt", "w");
	//ffff = fopen("vertex_stats.txt", "w");


	// loop pres vsechny castice:
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			rcon.compute_cell(c, j, i);		// spocte se bunka aktualne uvazovane castice

// CELL		(i,j; c)
			fprintf(f, "%d %g %g %g %d %d %g \n",			
				rcon.id[j][i],				// ID
											// bunka:
				c.volume(),					// objem
				c.surface_area(),			// povrch
				c.total_edge_distance(),	// delka hran
				c.number_of_faces(),		// pocet sten = pocet sousedu
				c.number_of_edges(),		// pocet hran
											// pocet vrcholu bunky = pocet hran - pocet sten + 2
				c.max_radius_squared()		// druha mocnina maximalni vzdalenosti vrcholu od stredu
			);

			c.neighbors(neigh);																													// F,E,V
			c.face_areas(fare);																													// F,E
			c.face_orders(ford);																												// F
			c.face_perimeters(fper);																											// F
			v = c.volume();																														// F

			edge_lengths(lengths, c);																											// E
			dihedral_angles(angles, c);																											// E

			for (k = 0; k < neigh.size(); k++) {	// loop over neighbors = loop over faces
				if (rcon.id[j][i] < neigh[k]) {		// if ID of neighbor is higher than id of current cell continue (to avoid duble counting)
					find_pos(ijk, q, neigh[k], &rcon);
					rcon.compute_cell(d, ijk, q);

// FACE		(i,j,k; c,d; ijk,q; neigh,ford,fare,fper; v;)
					fprintf(ff, "%g %d %g %g %g\n",
						fare[k],		// povrch steny
						ford[k],		// pocet vrcholu steny = pocet hran steny
						fper[k],		// obvod steny
						v,				// objem prvni bunky
						d.volume()		// objem druhe bunky
					);

					d.neighbors(neigh1);				// E
					d.face_areas(fare1);				// E
					dihedral_angles(angles1, d);		// E

					for (l = 0; l < neigh.size(); l++) {	// loop over neighbors again = loop over faces again
						if (neigh[k] < neigh[l]) {					// the second IF (to avoid triple counting) - every edge will appear only once 
							find_pos(ijk, q, neigh[l], &rcon);
							if (are_neighbors(d, ijk, q, &rcon)) {	// 3 neighboring cells determine edge; 

								// using the 1st cell:
								p = common_edge(c, k, l);			// find common edge in the 1st cell for faces of the 2nd and the 3rd cell, and determine its order p

// EDGE		(i,j,k,l; c,d,e; ijk,q; neigh,fare,neigh1,fare1,neigh2,lengths,angles,angles1,angles2; p,i1,i2)
								fprintf(fff, "%g %g %g ",
									lengths[p],		// delka hrany
									fare[k],		// obsah prvni steny, kt lze najit v prvni bunce
									fare[l]			// obsah druhe steny, kt lze najit v prvni bunce
									// dihedralni uhel v prvni bunce jiz umime urcit, ale nejdriv spocteme obsah treti steny
								);

								// using the 2nd cell:
								i1 = find_index(neigh1, rcon.id[j][i]);		// find the 1st and 3rd cell between neighbors of the 2nd cell
								i2 = find_index(neigh1, neigh[l]);

								fprintf(fff, "%g %g ", 
									fare1[i2],		// obsah treti steny
									angles[p]		// dihedralni uhel v radianech prislusejici prvni bunce c
								); 
								
								p = common_edge(d, i1, i2);

								fprintf(fff, "%g ", angles1[p]);  // dihedralni uhel prislusejici druhe bunce d

								// using the 3rd cell:
								rcon.compute_cell(e, ijk, q);		// E
								e.neighbors(neigh2);				// E
								dihedral_angles(angles2, e);		// E

								i1 = find_index(neigh2, rcon.id[j][i]);
								i2 = find_index(neigh2, neigh[k]);

								p = common_edge(e, i1, i2);

								fprintf(fff, "%g \n", angles2[p]);  // dihedralni uhel prislusejici treti bunce e

							}
						}
					}
				}
			}

// VERTEX - fce vertex stats vyuziva jiny pristup; dalo by se zde pokracovat, ale tezko rict zda/nakolik by to bylo vyhodne
			
		}
	}
	fclose(f);
	fclose(ff);
	fclose(fff);
}


// fc face_edge_lengths creates vector of lengths of edges in the same manner as the vector created by face_vertices
void face_edge_lengths(std::vector<double> &rv, voro::voronoicell_neighbor &rc)
{
	// [out]	v				the vector containing numbers of edges and their lenghts for every face of the cell in the same order as
	//							vector of vertices (vert)
	// [in]		c				the cell under consideration
	
	// observation: number of edges = number of vertices; for every face
	int i,j,k,l, v1,v2;
	double m;
	std::vector<int> rvert;

	rc.face_vertices(rvert);
	l = 0;

	rv.clear();
	rv.resize(rvert.size());

	for (i = 0; i < rc.number_of_faces(); i++) {		// loop over faces of the given cell c
		k = rvert[l];
		rv[l] = k;															// number of edges=vertices of this face
		for (j = 0; j < k-1; j++) {
			v1 = rvert[l + 1 + j]; v2 = rvert[l + 1 + j + 1];				// two vertices creating one edge
			
			m = sqrt( pow((rc.pts[3*v2] - rc.pts[3*v1]) * 0.5,2) + pow((rc.pts[3*v2+1] - rc.pts[3*v1+1]) * 0.5,2) + 
						pow((rc.pts[3*v2+2] - rc.pts[3*v1+2]) * 0.5,2));	// distance of these two vertices

			rv[l + 1 + j] = m;
		}

		v1 = rvert[l + k]; v2 = rvert[l + 1];
		m = sqrt(pow((rc.pts[3 * v2] - rc.pts[3 * v1]) * 0.5, 2) + pow((rc.pts[3 * v2 + 1] - rc.pts[3 * v1 + 1]) * 0.5, 2) + pow((rc.pts[3 * v2 + 2] - rc.pts[3 * v1 + 2]) * 0.5, 2));
		rv[l + k] = m;

		l = l + k + 1;
	}
}

// fc edge_lengths creates vector of lengths of edges ordered according to structure ed if we take into consideration only edges v1 < v2
void edge_lengths(std::vector<double> &rv, voro::voronoicell_neighbor &rc)
{
	// [out]	v				the vector of the lenghts of the edges in the specific order
	// [in]		c				the cell under consideration

	int i, j, k;
	k = 0;

	rv.clear();
	rv.resize(rc.number_of_edges());
	  
	for (i = 0; i < rc.p; i++) {			// loop over all vertices of the cell
		for (j = 0; j < rc.nu[i]; j++) {	// loop over all neighboring vertices of the given vertex
			if (i < rc.ed[i][j]) {			// skip neighboring vertices with the smaller number wrt ordering - to avoid double counting
				rv[k] = sqrt(pow((rc.pts[3*rc.ed[i][j]] - rc.pts[3*i]) * 0.5, 2) + pow((rc.pts[3*rc.ed[i][j] + 1] - rc.pts[3*i + 1]) * 0.5, 2) +
							pow((rc.pts[3 * rc.ed[i][j] + 2] - rc.pts[3 * i + 2]) * 0.5, 2));
							
				k++;
				
			}
		}
	}
}

void dihedral_angles(std::vector<double> &rv, voro::voronoicell_neighbor &rc)
{
	// [out]	v				the vector of the lenghts of the edges in the specific order
	// [in]		c				the cell under consideration

	int i, j, k, l, p, f1, f2, m,n,e;
	double a1, a2, a3, b1, b2, b3, c1, c2, c3, a,b,c,d, aa,bb,cc,dd, cosfi;
	std::vector<int> vert, ord;
	m = 0; n = 0; e = 0; 

	rv.clear();
	//rv.resize(rc.number_of_edges());

	rc.face_vertices(vert);

	for (i = 0; i < rc.p; i++) {			// loop over all vertices of the cell
		for (j = 0; j < rc.nu[i]; j++) {	// loop over all neighboring vertices of the given vertex
			if (i < rc.ed[i][j]) {			// skip neighboring vertices with the smaller number wrt ordering - to avoid double counting

			// najdi steny obsahujici tuto hranu:
				p = 0;
				for (k = 0; k < rc.number_of_faces(); k++) {
					for (l = 0; l < vert[p]; l++) {
						if (vert[p + 1 + l] == i) { m = 1; }				// je tam prvni vrchol
						if (vert[p + 1 + l] == rc.ed[i][j]) { n = 1; }		// je tam druhy vrchol
						
					}
					p = p + vert[p] + 1;
					if (m*n == 1) { if (e == 0) { f1 = k; e++; } else { f2 = k; } } // jsou tam oba
					m = 0; n = 0;
				}
				e = 0;
			// steny nalezeny

			// najdi rovnice rovin obsahujici tyto steny, tj nejdi tri body (vrcholy) z kazde steny:
				rc.face_orders(ord);
				// prvni strana:
				p = 0;
				for (k = 0; k < f1; k++) {
					p = p + ord[k];
				}
				p = p + f1 + 1;   // p nyni ukazuje na prvni vrchol f1-te steny
				// tri body:
				a1 = rc.pts[3 * vert[p]];
				a2 = rc.pts[3 * vert[p] + 1];
				a3 = rc.pts[3 * vert[p] + 2];
				b1 = rc.pts[3 * vert[p + 1]];
				b2 = rc.pts[3 * vert[p + 1] + 1];
				b3 = rc.pts[3 * vert[p + 1] + 2];
				c1 = rc.pts[3 * vert[p + 2]];
				c2 = rc.pts[3 * vert[p + 2] + 1];
				c3 = rc.pts[3 * vert[p + 2] + 2];
				// smerove vektory:
				b1 = b1 - a1;
				b2 = b2 - a2;
				b3 = b3 - a3;
				c1 = c1 - a1;
				c2 = c2 - a2;
				c3 = c3 - a3;
				// vektorovy soucin:
				a = b2*c3 - b3*c2;
				b = b3*c1 - b1*c3;
				c = b1*c2 - b2*c1;
				d = -a*a1 - b*a2 - c*a3;
				// druha strana:
				p = 0;
				for (k = 0; k < f2; k++) {
					p = p + ord[k];
				}
				p = p + f2 + 1;
				// tri body:
				a1 = rc.pts[3 * vert[p]];
				a2 = rc.pts[3 * vert[p] + 1];
				a3 = rc.pts[3 * vert[p] + 2];
				b1 = rc.pts[3 * vert[p + 1]];
				b2 = rc.pts[3 * vert[p + 1] + 1];
				b3 = rc.pts[3 * vert[p + 1] + 2];
				c1 = rc.pts[3 * vert[p + 2]];
				c2 = rc.pts[3 * vert[p + 2] + 1];
				c3 = rc.pts[3 * vert[p + 2] + 2];
				// smerove vektory:
				b1 = b1 - a1;
				b2 = b2 - a2;
				b3 = b3 - a3;
				c1 = c1 - a1;
				c2 = c2 - a2;
				c3 = c3 - a3;
				// vektorovy soucin:
				aa = b2*c3 - b3*c2;
				bb = b3*c1 - b1*c3;
				cc = b1*c2 - b2*c1;
				dd = -aa*a1 - bb*a2 - cc*a3;
			// mame rovnice rovin - a,b,c,d , aa,bb,cc,dd

			// spocteme dihedralni uhel techto rovin:
				cosfi = (a*aa + b*bb + c*cc) / (sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2))*sqrt(pow(aa, 2) + pow(bb, 2) + pow(cc, 2)));
				rv.push_back(acos(cosfi));
				//rv[o] = acos(cosfi);
				//o++;
				// acos - Returns the principal value of the arc cosine of x, expressed in radians.
			}
		}
	}

}

