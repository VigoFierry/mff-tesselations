
#include "Header.h"


using namespace voro; 




// fce feasibility overuje tzv. hardcore interactions; pro rozliseni nekolika specifickych situaci uvazujme pretizeni teto fce 

// fce feasibility pro pocatecni konfiguraci; verze pro container a container_poly
bool feasibility(container &con, const double &alfa, const double &beta, const double &B) 
{
	// [in]		con				the container with stored particles.
	// [in]	    alfa, beta		restrictions of the cell distances.
	// [in]	    B				restriction on the cell volume/shape.

	unsigned int k, fng;
	int ijk,q, i, j;
	int id;
	double x, y, z, xn, yn, zn, vol, voln, dist, dist2;
	voronoicell_neighbor c;  // bunka
	voronoicell d;
	std::vector<int> neigh;  // vektor pro prirazeni ID sousedu
	std::vector<int> vert;	 // vektor obsahujici vrcholy kazde steny bunky v nejakem danem poradi

	// ted chci projit pres vsechny castice, urcit jejich sousedy v mozaice, a pro ty s vyssim ID nez je samotne
	// ID cislo aktualni castice spocitat jejich vzdalenost od aktualni castice; pokud by tanto vzdalenost nebyla v danych
	// mezich, chci ukoncit cyklus a vratit FALSE; pokud cyklus projde bez ukonceni vsechny castice chci vratit TRUE
	// ---> tj chci vyuzit SYMETRIE

	// pro prochazeni lze pozit tridu loop z Voro++ - neni nezbytnosti, lze i jinak (cyklit pres boxy)
	// c_loop_all clo(con);          // loop pres vsechny castice
	// c_loop_order clo(con, po);    // loop jen pres vybrane

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			id = con.id[j][i];
			x = con.p[j][3 * i];
			y = con.p[j][3 * i + 1];
			z = con.p[j][3 * i + 2];

			con.compute_cell(c, j, i);

			if (B > 0) {
				vol = c.volume();	// compute its volume
			}
			c.neighbors(neigh);  // computes a vector list of neighbors of the current particle under consideration

			c.face_vertices(vert);		// computes list of vertices for each face of the cell
			fng = 0;

			// number of faces
			if (c.number_of_faces() > 12) { return false; };

			// loop over the neighbors
			for (k = 0; k < neigh.size(); k++) {

				// Skip if the neighbor information is smaller than this particle's ID, to avoid double counting. This also
				// removes faces that touch the walls, since the neighbor information is set to negative numbers for these cases. (neplati pri periodicite)
				if (neigh[k] > id) {

					face_dist(fng, vert, x, y, z, xn, yn, zn, c);		// {PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
					xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;	// {PER} urcim skutecne souradnice tohoto souseda za stenou

																		// feasibility musi kalkulovat se "skutecnymi" souradnicemi !!!!!!!!!!!!!
																		// stacilo by porovnat vzdalenosti steny od generatoru/stredu 
					dist = point_dist(x, y, z, xn, yn, zn);
					dist = dist / 2;

					//if (dist > beta) { std::cout << "Tesselation is NOT feasible (upper bound - " << dist << ").\n"; return false; }; 
					//if (dist < alfa) { std::cout << "Tesselation is NOT feasible (lower bound - " << dist << ").\n"; return false; };

					if (dist > beta) { return false; };
					if (dist < alfa) { return false; };


					if (B > 0) {  // B=0 znamena bez omezeni
						find_pos(ijk, q, neigh[k], &con);			 // najdu pozici aktualniho souseda 
						con.compute_cell(d, ijk, q);				 // spocte sousedni bunku
						voln = d.volume();							 // a jeji objem

						dist2 = pow(dist, 3);
						//if (dist2 > B*vol) { std::cout << "Tesselation is NOT feasible (shape constraint - " << (dist2 / vol) << " vs. " << B << ").\n"; return false; };
						//if (dist2 > B*voln) { std::cout << "Tesselation is NOT feasible (shape constraint - " << (dist2 / voln) << " vs. " << B << ").\n"; return false; };

						if (dist2 > B*vol) { return false; }; 
						if (dist2 > B*voln) { return false; };
					}
				}
				fng = fng + vert[fng] + 1;			// set actual position in vector of face vertices
			}
		}
	}
	return true;
}
/*
bool feasibility(container &con, const double &alfa, const double &beta, const double &B)
{
	// [in]		con				the container with stored particles.
	// [in]	    alfa, beta		restrictions of the cell distances.
	// [in]	    B				restriction on the cell volume/shape.

	// WARNING: vzdalenosti v POWER metrice nejsou symetricke !!!

	int i, j;
	int id;
	double x, y, z, vol, h_min, h_max;
	voronoicell_neighbor c;  // bunka
							 //std::vector<int> neigh;  // vektor pro prirazeni ID sousedu
	std::vector<int> vert;	 // vektor obsahujici vrcholy kazde steny bunky v nejakem danem poradi

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			id = con.id[j][i];
			x = con.p[j][4 * i];
			y = con.p[j][4 * i + 1];
			z = con.p[j][4 * i + 2];

			con.compute_cell(c, j, i);

			if (B > 0) {
				vol = c.volume();	// compute its volume
			}

			h_fcs(c, x, y, z, h_max, h_min);

			// alfa
			if (alfa > h_min) { std::cout << "Tesselation is NOT feasible (lower bound - " << h_min << ").\n"; return false; }
			//if (alfa > h_min) { return false; }

			// beta
			if (beta < h_max) { std::cout << "Tesselation is NOT feasible (upper bound - " << h_max << ").\n"; return false; }
			//if (beta < h_max) { return false; }

			// B
			if (B > 0) {	// B=0 means no constrain, i.e. in fact B=inf
				if ((pow(h_max, 3) / vol) > B) { std::cout << "Tesselation is NOT feasible (shape constraint - " << (pow(h_max, 3) / vol) << " vs. " << B << ").\n"; return false; };
				//if ((pow(h_max, 3)/vol) > B) { return false; };
			}

			// ...
		}
	}
	return true;
}*/

bool feasibility(container_poly &con, const double &alfa, const double &beta, const double &B)
{
	// [in]		con				the container with stored particles.
	// [in]	    alfa, beta		restrictions of the cell distances.
	// [in]	    B				restriction on the cell volume/shape.

	// WARNING: vzdalenosti v POWER metrice nejsou symetricke !!!

	int i, j;
	int id;
	double x, y, z, vol, h_min, h_max, tx, ty, tz;
	voronoicell_neighbor c;  // bunka
	bool grain;
	std::vector<int> vert;	 // vektor obsahujici vrcholy kazde steny bunky v nejakem danem poradi

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			id = con.id[j][i];
			x = con.p[j][4 * i];
			y = con.p[j][4 * i + 1];
			z = con.p[j][4 * i + 2];

			grain = con.compute_cell(c, j, i);

			if (grain == true) {	// the cell is nonempty

				c.centroid(tx, ty, tz);  // for Laguerre it makes a better sense to consider distance between faces and barycenter instead of generator

				tx = x + tx; ty = y + ty; tz = z + tz;		// real coordinates of barycenter

				if (B > 0) {
					vol = c.volume();	// compute its volume
				}

				h_fcs(c, x, y, z, tx, ty, tz, h_max, h_min);

				// alfa
				//if (alfa > h_min) { std::cout << "Tesselation is NOT feasible (lower bound - " << h_min << ").\n"; return false; }
				if (alfa > h_min) { return false; }

				// beta
				//if (beta < h_max) { std::cout << "Tesselation is NOT feasible (upper bound - " << h_max << ").\n"; return false; }
				if (beta < h_max) { return false; }

				// B
				if (B > 0) {	// B=0 means no constrain, i.e. in fact B=inf
					//if ((pow(h_max, 3) /  vol) > B) { std::cout << "Tesselation is NOT feasible (shape constraint - " << (pow(h_max, 3) / vol) << " vs. " << B << ").\n"; return false; };
					if ((pow(h_max, 3) / vol) > B) { return false; };
				}

				// ...

			} // END..if(nonempty cell)
		}
	}
	return true;
}


bool overlap_f(container_poly &con, const double &iota)
{
	// [in]		con				the container with stored particles.
	// [in]	    iota			restriction on the generators overlap .

	// WARNING: vzdalenosti v POWER metrice nejsou symetricke !!!

	double R_0 = 0.05; // (max. of radii)
	R_0 = R_0 + iota;

	int i, j, k, l;
	int id1, id2;
	double x, y, z, r, xx, yy, zz, rr, dist, a, b, c;
	bool xo = false, xoo = false, yo = false, yoo = false, zo = false, zoo = false;
	bool xxo = false, xxoo = false, yyo = false, yyoo = false, zzo = false, zzoo = false;
	int tou1, tou2;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			xo = false, xoo = false, yo = false, yoo = false, zo = false, zoo = false;

			id1 = con.id[j][i];
			x = con.p[j][4 * i];
			y = con.p[j][4 * i + 1];
			z = con.p[j][4 * i + 2];
			r = con.p[j][4 * i + 3];

			if (x - r < 0) { xo = true; }
			if (x + r > 1) { xoo = true; }
			if (y - r < 0) { yo = true; }
			if (y + r > 1) { yoo = true; }
			if (z - r < 0) { zo = true; }
			if (z + r > 1) { zoo = true; }

			for (l = 0; l < con.nxyz; l++) { // loop over boxes
				for (k = 0; k < con.co[l]; k++) { // loop over particles in considered box

					id2 = con.id[l][k];

					if (id1 < id2) {
						xxo = false, xxoo = false, yyo = false, yyoo = false, zzo = false, zzoo = false;
					
						xx = con.p[l][4 * k];
						yy = con.p[l][4 * k + 1];
						zz = con.p[l][4 * k + 2];
						rr = con.p[l][4 * k + 3];

						if (xx - rr < 0) { xxo = true; }
						if (xx + rr > 1) { xxoo = true; }
						if (yy - rr < 0) { yyo = true; }
						if (yy + rr > 1) { yyoo = true; }
						if (zz - rr < 0) { zzo = true; }
						if (zz + rr > 1) { zzoo = true; }

						// original image:
						dist = point_dist(x, y, z, xx, yy, zz);
						if ( dist < R_0) { // maximal distance R_0 + iota
						//	rr = con.p[l][4 * k + 3];

							if (iota < (dist-r)) { return false; }
							if (iota < (dist-rr)) { return false; }
						}					

						// periodic images:
						//		determine grain with higher number of touches to walls
						//		test all periodic images of this grain (max=8 in 3D)
						tou1 = xo + yo + zo + xoo + yoo + zoo;
						tou2 = xxo + yyo + zzo + xxoo + yyoo + zzoo;
						if (tou1 + tou2 > 0) {	// number of touches > 0
							if (tou1 > tou2) {	// the first grain has higher number of periodic images
								a = x; b = y; c = z;
								// ? possibilities:														JAK OPTIMALIZOVAT ???
								if (xo == true) { a = x - 1; }
								if (yo == true) { b = y - 1; }
								if (zo == true) { c = z - 1; }
								if (xoo == true) { a = x + 1; }
								if (yoo == true) { b = y + 1; }
								if (zoo == true) { c = z + 1; }

								// 7 possibilities:
								dist = point_dist(a, y, z, xx, yy, zz);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
								dist = point_dist(x, b, z, xx, yy, zz);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
								dist = point_dist(x, y, c, xx, yy, zz);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
								dist = point_dist(a, b, z, xx, yy, zz);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
								dist = point_dist(a, y, c, xx, yy, zz);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
								dist = point_dist(x, b, c, xx, yy, zz);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
								dist = point_dist(a, b, c, xx, yy, zz);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
								/*
								if (xo == true){ 
									a = x - 1;
									dist = point_dist(a, y, z, xx, yy, zz);
									if (iota < (dist - r)) { return false; }
									if (iota < (dist - rr)) { return false; }
								}
								if (yo == true) {
									b = y - 1;
									dist = point_dist(x, b, z, xx, yy, zz);
									if (iota < (dist - r)) { return false; }
									if (iota < (dist - rr)) { return false; }
								}
								if (zo == true) {
									c = z - 1;
									dist = point_dist(x, y, c, xx, yy, zz);
									if (iota < (dist - r)) { return false; }
									if (iota < (dist - rr)) { return false; }
								}
								if (xoo == true) {
									a = x + 1;
									dist = point_dist(a, y, z, xx, yy, zz);
									if (iota < (dist - r)) { return false; }
									if (iota < (dist - rr)) { return false; }
								}
								if (yoo == true) {
									b = y + 1;
									dist = point_dist(x, b, z, xx, yy, zz);
									if (iota < (dist - r)) { return false; }
									if (iota < (dist - rr)) { return false; }
								}
								if (zoo == true) {
									c = z + 1;
									dist = point_dist(x, y, c, xx, yy, zz);
									if (iota < (dist - r)) { return false; }
									if (iota < (dist - rr)) { return false; }
								}
								if (zo == true) {
									c = z - 1;
									dist = point_dist(x, y, z, xx, yy, zz);
									if (iota < (dist - r)) { return false; }
									if (iota < (dist - rr)) { return false; }
								}
								if (zo == true) {
									c = z - 1;
									dist = point_dist(x, y, z, xx, yy, zz);
									if (iota < (dist - r)) { return false; }
									if (iota < (dist - rr)) { return false; }
								}
								if (zo == true) {
									c = z - 1;
									dist = point_dist(x, y, z, xx, yy, zz);
									if (iota < (dist - r)) { return false; }
									if (iota < (dist - rr)) { return false; }
								}
								if (zo == true) {
									c = z - 1;
									dist = point_dist(x, y, z, xx, yy, zz);
									if (iota < (dist - r)) { return false; }
									if (iota < (dist - rr)) { return false; }
								}*/
							}
							else {				// the second grain has higher number of periodic images
								a = xx; b = yy; c = zz;
								// ? possibilities:
								if (xxo == true) { a = xx - 1; }
								if (yyo == true) { b = yy - 1; }
								if (zzo == true) { c = zz - 1; }
								if (xxoo == true) { a = xx + 1; }
								if (yyoo == true) { b = yy + 1; }
								if (zzoo == true) { c = zz + 1; }

								// 7 possibilities:
								dist = point_dist(x, y, z, a, yy, zz);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
								dist = point_dist(x, y, z, xx, b, zz);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
								dist = point_dist(x, y, z, xx, yy, c);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
								dist = point_dist(x, y, z, a, b, zz);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
								dist = point_dist(x, y, z, a, yy, c);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
								dist = point_dist(x, y, z, xx, b, c);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
								dist = point_dist(x, y, z, a, b, c);
								if (iota < (dist - r)) { return false; }
								if (iota < (dist - rr)) { return false; }
							}
							
						}
						
					} // end if(id)

				}
			} // end second loop		
		}
	} // end first loop
	return true;
}

// fce feasibility pro podokno, verze pouze pro container_poly; 
bool feasibility_subwindows(container_poly &con, const double &alfa, const double &beta, const double &B, double x, double y, double z, double rad)
{
	// [in]		con				the container with stored particles.
	// [in]	    alfa, beta		restrictions of the cell distances.
	// [in]	    B				restriction on the cell volume/shape.
	// [in]		x,y,z			coordinates of particle under examination
	// [in]		rad				its radius

	// WARNING: vzdalenosti v POWER metrice nejsou symetricke !!!

	int i, j;
	int id;
	double xx, yy, zz, vol, h_min, h_max, tx, ty, tz;
	voronoicell_neighbor c;  // bunka
	//std::vector<int> neigh;  // vektor pro prirazeni ID sousedu
	std::vector<int> vert;	 // vektor obsahujici vrcholy kazde steny bunky v nejakem danem poradi


	// translation: 
	double trans = 0.4;	// suitable choice of translation is very important (it could be based on the radius value)

	double lx, ux, ly, uy, lz, uz;
	int t1, t2, t3, t;
	bool test;

	lx = x - trans; ux = x + trans;
	ly = y - trans; uy = y + trans;
	lz = z - trans; uz = z + trans;

	t1 = 0; if (lx < 0) { t1 = -1; }
	else { if (ux > 1) { t1 = 1; } }
	t2 = 0; if (ly < 0) { t2 = -1; }
	else { if (uy > 1) { t2 = 1; } }
	t3 = 0; if (lz < 0) { t3 = -1; }
	else { if (uz > 1) { t3 = 1; } }

	t = t1 + t2 + t3;


	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			id = con.id[j][i];
			xx = con.p[j][4 * i];
			yy = con.p[j][4 * i + 1];
			zz = con.p[j][4 * i + 2];

			test = false; // indicator which parlicle examine or not

			// testing conditions:

			if (lx < xx && xx < ux && ly < yy && yy < uy && lz < zz && zz < uz) { test = true; }	// particle in the original subwindow
			// translations in one coordinate:
			if (t1 == 0) {} else { if (lx - t1 < xx && xx < ux - t1 && ly < yy && yy < uy && lz < zz && zz < uz) { test = true; }; }
			if (t2 == 0) {} else { if (lx < xx && xx < ux && ly - t2 < yy && yy < uy - t2 && lz < zz && zz < uz) { test = true; }; }
			if (t3 == 0) {} else { if (lx < xx && xx < ux && ly < yy && yy < uy && lz - t3 < zz && zz < uz - t3) { test = true; }; }
			// translations in two coordinates:
			if (t1*t2 == 0) {} else { if (lx - t1 < xx && xx < ux - t1 && ly - t2 < yy && yy < uy - t2 && lz < zz && zz < uz) { test = true; } }
			if (t1*t3 == 0) {} else { if (lx - t1 < xx && xx < ux - t1 && ly < yy && yy < uy && lz - t3 < zz && zz < uz - t3) { test = true; } }
			if (t2*t3 == 0) {} else { if (lx < xx && xx < ux && ly - t2 < yy && yy < uy - t2 && lz - t3 < zz && zz < uz - t3) { test = true; } }
			// translation in three coordinates:
			if (t1*t2*t3 == 0) {} else { if (lx - t1 < xx && xx < ux - t1 && ly - t2 < yy && yy < uy - t2 && lz - t3 < zz && zz < uz - t3) { test = true; } }

			if (test == true) {
				con.compute_cell(c, j, i);
				c.centroid(tx, ty, tz);

				tx = xx + tx; ty = yy + ty; tz = zz + tz;	// real coordinates of centroid

				if (B > 0) {
					vol = c.volume();	// compute its volume
				}

				h_fcs(c, xx, yy, zz, tx, ty, tz, h_max, h_min);

				// alfa
				//if (alfa > h_min) { std::cout << "Tesselation is NOT feasible (lower bound - " << h_min << ").\n"; return false; }
				if (alfa > h_min) { return false; }

				// beta
				//if (beta < h_max) { std::cout << "Tesselation is NOT feasible (upper bound - " << h_max << ").\n"; return false; }
				if (beta < h_max) { return false; }

				// B
				if (B > 0) {	// B=0 means no constrain, i.e. in fact B=inf
					//if ((pow(h_max, 3) / vol) > B) { std::cout << "Tesselation is NOT feasible (shape constraint - " << (pow(h_max, 3) / vol) << " vs. " << B << ").\n"; return false; };
					if ((pow(h_max, 3)/vol) > B) { return false; };
				}

				// ...

			} // END{if; test==true)
		}
	}
	return true;
}



// fce feasibility pro dany seznam bodu, verze pro container a container_poly 
bool feasibility(voro::container &con, std::vector<int> cells, const double &alfa, const double &beta, const double &B)
// argumentem ma byt seznam bodu, tj. array, predem nespecifikovane delky - vektor poloh castic v containeru, fci musime
//    navic predat i cely container, kde jsou ulozene pozice castic
{
	// [in]		con				the container with stored particles.
	// [in]		cells			the list of particles which should be verified.
	// [in]	    alfa, beta		restrictions of the cell distances.

	unsigned int i;
//	int n;
	double x, y, z, h_min, h_max, vol;
	voronoicell_neighbor c;  // bunka
	std::vector<int> vert;	 // vektor obsahujici vrcholy kazde steny bunky v nejakem danem poradi

	for (i = 0; i < cells.size() / 2; i++) {				// take and examine each cell separatedly
		con.compute_cell(c, cells[2 * i], cells[2 * i + 1]);

		x = con.p[cells[2 * i]][3*cells[2 * i + 1] ];
		y = con.p[cells[2 * i]][3*cells[2 * i + 1] + 1];
		z = con.p[cells[2 * i]][3*cells[2 * i + 1] + 2];
		
		vol = c.volume();  // potreba pro B
		
		h_fcs(c, x, y, z, x, y, z, h_max, h_min);

		// alfa
		//if (h_min < alfa) { std::cout << "Tesselation is NOT feasible (lower bound - " << h_min << ").\n"; return false; };
		if (h_min < alfa) { return false; };

		// beta
		//if (h_max > beta) { std::cout << "Tesselation is NOT feasible (upper bound - " << h_max << ").\n"; return false; }; 
		if (h_max > beta) { return false; }; 
			
		// B
		if (B > 0) {	// B=0 means no constrain, i.e. in fact B=inf
			// note: use the appropriate power with respect to dimension, i.e. power 2 in 2D, power 3 in 3D, ...
			//if ((pow(h_max, 3) / vol) > B) { std::cout << "Tesselation is NOT feasible (shape constraint - " << (pow(h_max, 3) / vol) << " vs. " << B << ").\n"; return false; };
			if ((pow(h_max, 3)/vol) > B) { return false; };
		}

		// number of faces
		//    if (c.number_of_faces() > 12) { return false; };

		// ...
	}
	return true;
}

bool feasibility(voro::container_poly &con, std::vector<int> cells, const double &alfa, const double &beta, const double &B)
// argumentem ma byt seznam bodu, tj. array, predem nespecifikovane delky - vektor poloh castic v containeru, fci musime
//    navic predat i cely container, kde jsou ulozene pozice castic
{
	// [in]		con				the container with stored particles.
	// [in]		cells			the list of particles which should be verified.
	// [in]	    alfa, beta		restrictions of the cell distances.

	// WARNING: power metric is not symetric

	unsigned int i;
//	int n;
	double x, y, z, vol, h_min, h_max, tx, ty, tz;
	bool grain;
	voronoicell_neighbor c;  // bunka
	std::vector<int> vert;	 // vektor obsahujici vrcholy kazde steny bunky v nejakem danem poradi

	for (i = 0; i < cells.size() / 2; i++) {				// take and examine each cell separatedly
		if (cells[2 * i] == -1) { grain = false; }			// je-li cells[2 * i] == -1, tak i cells[2 * i + 1] == -1, a znamena to, ze i-ta castice byla odstranena
		else { grain = con.compute_cell(c, cells[2 * i], cells[2 * i + 1]); }

		if (grain == true) {		// nonempty cell
			x = con.p[cells[2 * i]][4 * cells[2 * i + 1]];
			y = con.p[cells[2 * i]][4 * cells[2 * i + 1] + 1];
			z = con.p[cells[2 * i]][4 * cells[2 * i + 1] + 2];

			vol = c.volume();  // potreba pro B
			c.centroid(tx, ty, tz);

			tx = x + tx; ty = y + ty; tz = z + tz;	// real coordinates of centroid

			h_fcs(c, x, y, z, tx, ty, tz, h_max, h_min);

			// alfa
			//if (h_min < alfa) { std::cout << "Tesselation is NOT feasible (lower bound - " << h_min << ").\n"; return false; };
			if (h_min < alfa) { return false; };

			// beta
			//if (h_max > beta) { std::cout << "Tesselation is NOT feasible (upper bound - " << h_max << ").\n"; return false; };
			if (h_max > beta) { return false; };

			// B
			if (B > 0) {	// B=0 means no constrain, i.e. in fact B=inf
				//if ((pow(h_max, 3)/ vol) > B) { std::cout << "Tesselation is NOT feasible (shape constraint - " << (pow(h_max, 3)/ vol) << " vs. " << B << ").\n"; return false; };
				if ((pow(h_max, 3) / vol) > B) { return false; };
			}

			// ...

		} // END..if(nonempty cell)
	}

	return true;
}


// fce feasibility pro specialni pripad pridani bodu 
bool feasibility(voro::container &con, int ijk1, int q1, const double &alfa, const double &beta, const double &B)
{
	// [in]		con				the container with stored particles.
	// [in]		ijk1, q1		the positions of the new added particle.
	// [in]	    alfa, beta		restrictions of the cell distances.
	// [in]	    B				restriction on the cell volume/shape.

	unsigned int i, j, fng;
	int ijk, q;
	double x, y, z, xn, yn, zn, vol, h_max2, dist;
	voronoicell_neighbor c,d;  // bunka s informaci o sousedech
	std::vector<int> neigh,neighi,vert, verti;  // vektor pro prirazeni ID sousedu

	// fce dostala primo udaje o castici (spec pripad BIRTH - jedna nova castice)

	con.compute_cell(c, ijk1, q1);     // Computes cell for given particle
	x = con.p[ijk1][3 * q1]; 
	y = con.p[ijk1][3 * q1 + 1]; 
	z = con.p[ijk1][3 * q1 + 2];       // Computes coordinates of this particle
	
    c.neighbors(neigh);              // Computes a vector list of neighbors of the current particle under consideration
	c.face_vertices(vert);			 // Computes list of vertices for each face of the cell

	if (B > 0) {
		vol = c.volume();	// compute its volume
	}

	fng = 0;

							 // loop over the neighbors
	for (i = 0; i<neigh.size(); i++) {

		face_dist(fng, vert, x, y, z, xn, yn, zn, c);		// {PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
		xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;	// {PER} urcim skutecne souradnice tohoto souseda za stenou
															// feasibility musi kalkulovat se "skutecnymi" souradnicemi !!!!!!!!!!!!!
		dist = point_dist(x, y, z, xn, yn, zn) / 2;

		if (dist > beta) { /*std::cout << "Tesselation is NOT feasible (upper bound - " << dist << ").\n";*/ return false; }; 
		if (dist < alfa) { /*std::cout << "Tesselation is NOT feasible (lower bound - " << dist << ").\n";*/ return false; };
	
		if (B > 0) {  // B=0 znamena bez omezeni

			if (pow(dist, 3) > B*vol) { /*std::cout << "Tesselation is NOT feasible (shape constraint - " << pow(dist, 3)/vol << " vs. " << B << ").\n";*/ return false; }; // double sqrt(double x);
			
			// pridanim bodu se zmensi objemy jeho sousedu, je proto nutne overit podminku pro vsechny bunky, jejichz objem se zmenil:
			find_pos(ijk, q, neigh[i], &con);			 // najdu pozici aktualniho souseda 
			con.compute_cell(d, ijk, q);				 // spocte jeho bunku
			d.neighbors(neighi);						 // a jeho sousedy
			d.face_vertices(verti);						 // computes list of vertices for each face of the cell
			vol = d.volume();							 // a jeji objem
			j = neighi.size();
			h_max2 = pow(h_maximum(j, verti, d, xn, yn, zn), 3);
			
			if (h_max2 > B*vol) { /*std::cout << "Tesselation is NOT feasible (shape constraint - " << h_max2 / vol << " vs. " << B << ").\n";*/ return false; }; // double sqrt(double x);
		}

		fng = fng + vert[fng] + 1;			// set actual position in vector of face vertices
	
	}

	return true;
}


// fce feasibility pro specialni pripad smazani bodu 
bool feasibility(voro::container &con, std::vector<int> &rsr, std::vector<int> &rsio, std::vector<double> &rsap, const double &alfa, const double &beta, const double &B)
{
	// [in]		con				the container with stored particles.
	// [in]		ijk1, q1		the positions of the new added particle.
	// [in]	    alfa, beta		restrictions of the cell distances.
	// [in]	    B				restriction on the cell volume/shape.

	unsigned int i, j;
	int ijk_i, q_i, ijk_j, q_j;
	double x, y, z, xn, yn, zn, vol, voln, dist, dist2;
	voronoicell_neighbor c;  // bunka s informaci o sousedech
	voronoicell d;

	for (i = 0; i < rsr.size() / 2; i++) {                           // {F} loop over secondary particles
		ijk_i = rsr[2 * i]; q_i = rsr[2 * i + 1];
		con.compute_cell(c, ijk_i, q_i);							// {F} compute cell of i-th secondary particle

		if (rsio[i] == 0) {
			x = con.p[ijk_i][3 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][3 * q_i + 1];
			z = con.p[ijk_i][3 * q_i + 2];
		}
		else {
			x = rsap[3 * (rsio[i] - 1)];
			y = rsap[3 * (rsio[i] - 1) + 1];
			z = rsap[3 * (rsio[i] - 1) + 2];
		}

		for (j = i + 1; j < rsr.size() / 2; j++) {                   // {F} loop over secondary particles with "higher order" (to avoid double counting) 
			ijk_j = rsr[2 * j]; q_j = rsr[2 * j + 1];
			if (are_neighbors(c, ijk_j, q_j, &con)) {	// {F} sousedi-li bunky i a j ... (i.e. if c_i ~ c_j)
				
				// std::cout << "are \n"; ///////////////////////////////////////////////////////////////////////////////
				if (rsio[j] == 0) {
					xn = con.p[ijk_j][3 * q_j];				// {F} computes coordinates of j-th neighbor
					yn = con.p[ijk_j][3 * q_j + 1];
					zn = con.p[ijk_j][3 * q_j + 2];
				}
				else {
					xn = rsap[3 * (rsio[j] - 1)];
					yn = rsap[3 * (rsio[j] - 1) + 1];
					zn = rsap[3 * (rsio[j] - 1) + 2];
				}
				// {F} vzdalenost dvou bodu: sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
				dist = point_dist(x, y, z, xn, yn, zn) / 2;

				if (dist > beta) { /*std::cout << "Tesselation is NOT feasible (" << dist << " - upper bound).\n";*/ return false; }; // double sqrt(double x);
				if (dist < alfa) { /*std::cout << "Tesselation is NOT feasible (" << dist << " - lower bound).\n";*/ return false; };
				//if (point_dist(x, y, z, xn, yn, zn) < alfa) { std::cout << "D: LOW.\n"; };

				if (B > 0) {  // B=0 znamena bez omezeni
					// odebranim bodu se zvetsi objem jeho sousedu, proto staci overit podminku jen pro nove vznikle relace sousednosti
					con.compute_cell(d, ijk_j, q_j);				 // spocte sousedni bunku
					vol = c.volume();
					voln = d.volume();							 // a jeji objem

					dist2 = pow(dist, 3);
					if (dist2 > B*vol) { /*std::cout << "Tesselation is NOT feasible (shape constraint - " << dist2 / vol << " vs. " << B << ").\n";*/ return false; }; // double sqrt(double x);
					if (dist2 > B*voln) { /*std::cout << "Tesselation is NOT feasible (shape constraint - " << dist2 / voln << " vs. " << B << ").\n";*/ return false; };
				}
			}
		}
	}

	return true;
}
// pozn.: fce feasibility pro add a delete nelze adaptovat pro pripad poly - kvuli nesymetricite nedojde ke zjednoduseni fce fesibility pro seznam bodu

/*
// fce BDMA_step - telo algoritmu, provadi operace ADD, DELETE, MOVE a vola prislusne funkce TRY pro spocteni prislusnych psti
void bdma_step(long &npart, std::vector<int> &fid, voro::container &con, double sigma, double theta, double zet, const double &alfa, const double &beta, const double &B, long &rn_add, long &rn_del, long &rn_mov)
{
	// [in,out]		npart	number of added particles in the whole history.
	// [in,out]		fid		vector of available id큦.
	// [in]			con		the container with stored particles.
	// [in]			sigma	parameter of normal distribution.
	// [in]			theta	parameter of V2 function
	// [in]			zet     intensity constant of reference Poisson process
	
	double rn1 = uniform(0, 1);   // random number between 0 and 1 
	double rn2 = uniform(0,1);   // random number between 0 and 1 
	//rn1 = 0.8;
	//rn2 = 0.0000000000000001;
	double nx, ny, nz, x, y, z;
	double pst; 
	double const1 = 1;
	double const2 = 2;
	double const3 = 3;
	
	int del, id, no;
	int ijk, q; // ijk_add, q_add; unsigned int i,j,k; 

	//double alfa, beta; // asi je treba predat teto funkci  // nebo je mit jako globalni promenne
	//alfa = 0.02;
	//beta = 0.2;

	//std::cout << "pst: " << rn1 << " " << rn2 << "\n";
	const1 = const1 / const3;
	const2 = const2 / const3;

	//std::cout << "fid: ";	///////////////////////////////////////////////////////////////////////////////////////////
	//for (i = 0; i < fid.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	//std::cout << fid[i] << " ";
	//} std::cout << " \n";

	if (rn1 <= (const1)) {
		nx = uniform(0, 1); ny = uniform(0, 1); nz = uniform(0, 1); // coordinates of new particle
		// what about ID ??? - dve moznosti - pamatovat pocet jiz pridanych castic, tj. pamatovat nejvyssi jiz pouzite ID
	    //                                  - vezt zaznam o uvolnenych ID po odstranenych casticich a znovu je pouzivat
		if (fid.size() > 0) {				// najiti vhodneho ID
			id = fid[0]; 			
		}   // zamerem bylo zjistit zdali je ve fid ulozene nejake volne id ci ne???
		else { id = npart + 1;}  // pokud nebylo dostupne ID ve fid, pouzije se maximalni hodnota ID + 1
		std::cout << "ADD " << id << " " << nx << " " << ny << " " << nz << " \n"; // popis cinosti

		no = con.total_particles();                             // pocet castic v kontejneru
		pst = try_add(id, nx, ny, nz, con, theta, alfa, beta, B);		  	// spocte pravdepodobnost se kterou dojde k operaci ADD
		// pst = 0.5;
		std::cout << pst << " ";  
		pst = pst*(zet / (no + 1));   
		std::cout << pst << " ";    

		// rn = rnd();							// generate random number between 0 and 1 
		if (rn2 < pst) {
			con.put(id, nx, ny, nz);		// pridani castice
			std::cout << "YES \n";
			if (id > npart) {				// pouzilo-li se ID o jedna vyssi nez maximalni ID v kontejneru a pridani
				npart = npart + 1;			// castice se uskutecni, pak je nutne navysit toto maximalni ID o 1
			}
			else {							// pouzilo-li se ID z fid, musi ve fid smazat
				fid.erase(fid.begin());
			}		
			rn_add++;						// counter of added particles
		}
		else { std::cout << "NO \n"; }
		//std::cout << npart << "\n";  /////////////////////////////////////////////////////////////////////////////////
	}     // BIRTH
	
	if (((const1) < rn1) && (rn1 <= (const2))) {
		del = uniform_int(1, con.total_particles());	// vybere castici (resp. jeji poradi v containeru) ke smazani (uniforme z celkoveho poctu castic)
		// del neni ID!!!
		find_part(ijk, q, del, &con);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = con.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		std::cout << "DELETE " << del << " " << con.p[ijk][3 * q] << " " << con.p[ijk][3 * q + 1] << " " << con.p[ijk][3 * q + 2] << "\n"; // popis cinosti
		// nyni je potreba tuto castici odstranit z datovych struktur containeru, nestane se container nestabilni???
		// pred vymazanim je treba si zapsat jeji sousedy:

		no = con.total_particles();                     // pocet castic v kontejneru
		pst = try_delete(ijk, q, con, theta, alfa, beta, B);		// spocte pravdepodobnost se kterou dojde k operaci DELETE
		// pst = 0.5;
		std::cout << pst << " ";
		pst = pst*(no / zet);
		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {
			erase(ijk, q, &fid, &con);						// smazani castice - uvolni ID do fid k opetovnemu pouziti
			std::cout << "YES \n";
			rn_del++;										// counter of deleted particles
		}
		else { std::cout << "NO \n"; }
	}     // DEATH
	
	if (rn1 > (const2)) {
		del = uniform_int(1, con.total_particles());	// vybere castici ke smazani (uniforme z celkoveho poctu castic)
		// del neni ID!!!
		find_part(ijk, q, del, &con);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = con.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = con.p[ijk][3 * q]; y = con.p[ijk][3 * q + 1]; z = con.p[ijk][3 * q + 2]; // urci jeji souradnice
		// vygeneruje novou castici z normalniho rozdeleni:
		nx = normal(x, sigma); ny = normal(y, sigma); nz = normal(z, sigma);   // coordinates of new particle
		//std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
		// these coordinates could be outside the window [0,1]^3 and in this case it is necessary to remap it into window:
		nx = nx - step_int(nx); ny = ny - step_int(ny); nz = nz - step_int(nz);
		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti

		pst = try_move(ijk, q, nx, ny, nz, con, theta, alfa, beta, B);		// urci pravdepodobnost se kterou dojde k operaci MOVE
		// pst = 0.5;
		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {
			erase(ijk, q, &con);					// smazani castice
			con.put(del, nx, ny, nz);				// pridani nove castice (ID zustava zachovano)
			 
			std::cout << "YES \n";
			rn_mov++;								// counter of moved particles
		}
		else { std::cout << "NO \n"; }
	}     // MOVE 
}*/


// fce BDMA_step - bez vypisovani
void bdma_step(long &npart, std::vector<int> &fid, voro::container &con, double sigma, double theta, double zet, const double &alfa, const double &beta, const double &B, long &rn_add, long &rn_del, long &rn_mov)
{
	// [in,out]		npart	number of added particles in the whole history.
	// [in,out]		fid		vector of available id큦.
	// [in]			con		the container with stored particles.
	// [in]			sigma	parameter of normal distribution.
	// [in]			theta	parameter of V2 function
	// [in]			zet     intensity constant of reference Poisson process

	double rn1 = uniform(0, 1);   // random number between 0 and 1 
	double rn2 = uniform(0, 1);   // random number between 0 and 1 
								  //rn1 = 0.8;
								  //rn2 = 0.0000000000000001;
	double nx, ny, nz, x, y, z;
	double pst;
	double const1 = 1;
	double const2 = 2;
	double const3 = 3;

	int del, id, no;
	int ijk, q; // ijk_add, q_add; unsigned int i,j,k; 

				//double alfa, beta; // asi je treba predat teto funkci  // nebo je mit jako globalni promenne
				//alfa = 0.02;
				//beta = 0.2;

				//std::cout << "pst: " << rn1 << " " << rn2 << "\n";
	const1 = const1 / const3;
	const2 = const2 / const3;

	if (rn1 <= (const1)) {
		nx = uniform(0, 1); ny = uniform(0, 1); nz = uniform(0, 1); // coordinates of new particle
																	// what about ID ??? - dve moznosti - pamatovat pocet jiz pridanych castic, tj. pamatovat nejvyssi jiz pouzite ID
																	//                                  - vezt zaznam o uvolnenych ID po odstranenych casticich a znovu je pouzivat
		if (fid.size() > 0) {				// najiti vhodneho ID
			id = fid[0];
		}   // zamerem bylo zjistit zdali je ve fid ulozene nejake volne id ci ne???
		else { id = npart + 1; }  // pokud nebylo dostupne ID ve fid, pouzije se maximalni hodnota ID + 1
//		std::cout << "ADD " << id << " " << nx << " " << ny << " " << nz << " \n"; // popis cinosti

		no = con.total_particles();                             // pocet castic v kontejneru
		pst = try_add(id, nx, ny, nz, con, theta, alfa, beta, B);		  	// spocte pravdepodobnost se kterou dojde k operaci ADD
																			// pst = 0.5;
//		std::cout << pst << " ";
		pst = pst*(zet / (no + 1));
//		std::cout << pst << " ";

		// rn = rnd();							// generate random number between 0 and 1 
		if (rn2 < pst) {
			con.put(id, nx, ny, nz);		// pridani castice
//			std::cout << "YES \n";
			if (id > npart) {				// pouzilo-li se ID o jedna vyssi nez maximalni ID v kontejneru a pridani
				npart = npart + 1;			// castice se uskutecni, pak je nutne navysit toto maximalni ID o 1
			}
			else {							// pouzilo-li se ID z fid, musi ve fid smazat
				fid.erase(fid.begin());
			}
			rn_add++;						// counter of added particles
		}
//		else { std::cout << "NO \n"; }
		//std::cout << npart << "\n";  /////////////////////////////////////////////////////////////////////////////////
	}     // BIRTH

	if (((const1) < rn1) && (rn1 <= (const2))) {
		del = uniform_int(1, con.total_particles());	// vybere castici (resp. jeji poradi v containeru) ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &con);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = con.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
//		std::cout << "DELETE " << del << " " << con.p[ijk][3 * q] << " " << con.p[ijk][3 * q + 1] << " " << con.p[ijk][3 * q + 2] << "\n"; // popis cinosti
																																		   // nyni je potreba tuto castici odstranit z datovych struktur containeru, nestane se container nestabilni???
																																		   // pred vymazanim je treba si zapsat jeji sousedy:

		no = con.total_particles();                     // pocet castic v kontejneru
		pst = try_delete(ijk, q, con, theta, alfa, beta, B);		// spocte pravdepodobnost se kterou dojde k operaci DELETE
																	// pst = 0.5;
//		std::cout << pst << " ";
		pst = pst*(no / zet);
//		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {
			erase(ijk, q, &fid, &con);						// smazani castice - uvolni ID do fid k opetovnemu pouziti
//			std::cout << "YES \n";
			rn_del++;										// counter of deleted particles
		}
//		else { std::cout << "NO \n"; }
	}     // DEATH

	if (rn1 > (const2)) {
		del = uniform_int(1, con.total_particles());	// vybere castici ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &con);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = con.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = con.p[ijk][3 * q]; y = con.p[ijk][3 * q + 1]; z = con.p[ijk][3 * q + 2]; // urci jeji souradnice
																					 // vygeneruje novou castici z normalniho rozdeleni:
		nx = normal(x, sigma); ny = normal(y, sigma); nz = normal(z, sigma);   // coordinates of new particle
																			   //std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
																			   // these coordinates could be outside the window [0,1]^3 and in this case it is necessary to remap it into window:
		nx = nx - step_int(nx); ny = ny - step_int(ny); nz = nz - step_int(nz);
//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti

		pst = try_move(ijk, q, nx, ny, nz, con, theta, alfa, beta, B);		// urci pravdepodobnost se kterou dojde k operaci MOVE
																			// pst = 0.5;
//		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {
			erase(ijk, q, &con);					// smazani castice
			con.put(del, nx, ny, nz);				// pridani nove castice (ID zustava zachovano)

//			std::cout << "YES \n";
			rn_mov++;								// counter of moved particles
		}
//		else { std::cout << "NO \n"; }
	}     // MOVE 
}

/*
void bdma_step_2(long &npart, std::vector<int> &fid, voro::container &con, double sigma, double theta, double zet, const double &alfa, const double &beta, const double &B, long &rn_add, long &rn_del, long &rn_mov, bool &fc)
{
	// [in,out]		npart	number of added particles in the whole history.
	// [in,out]		fid		vector of available id큦.
	// [in]			con		the container with stored particles.
	// [in]			sigma	parameter of normal distribution.
	// [in]			theta	parameter of V2 function
	// [in]			zet     intensity constant of reference Poisson process
	// [in]			fc		identificator of energy (0 - uses face areas, 1 - uses dihedral angles)

	double rn1 = uniform(0, 1);   // random number between 0 and 1 
	double rn2 = uniform(0, 1);   // random number between 0 and 1 
	//rn1 = 0.8; 
	//rn2 = 0.0000000000000001;
	double nx, ny, nz, x, y, z;
	double pst;
	double const1 = 1;
	double const2 = 2;
	double const3 = 3;

	int del, id, no;
	int ijk, q; // ijk_add, q_add; unsigned int i,j,k; 


	//std::cout << "pst: " << rn1 << " " << rn2 << "\n";
	const1 = const1 / const3;
	const2 = const2 / const3;

	if (rn1 <= (const1)) {
		nx = uniform(0, 1); ny = uniform(0, 1); nz = uniform(0, 1); // coordinates of new particle
																	// what about ID ??? - dve moznosti - pamatovat pocet jiz pridanych castic, tj. pamatovat nejvyssi jiz pouzite ID
																	//                                  - vezt zaznam o uvolnenych ID po odstranenych casticich a znovu je pouzivat
		if (fid.size() > 0) {				// najiti vhodneho ID
			id = fid[0];
		}   // zamerem bylo zjistit zdali je ve fid ulozene nejake volne id ci ne???
		else { id = npart + 1; }  // pokud nebylo dostupne ID ve fid, pouzije se maximalni hodnota ID + 1
		std::cout << "ADD " << id << " " << nx << " " << ny << " " << nz << " \n"; // popis cinosti

		no = con.total_particles();                             // pocet castic v kontejneru
		pst = try_add_2(id, nx, ny, nz, con, theta, alfa, beta, B, fc);		  	// spocte pravdepodobnost se kterou dojde k operaci ADD
																				// pst = 0.5;
		std::cout << pst << " ";
		pst = pst*(zet / (no + 1));
		std::cout << pst << " ";

		// rn = rnd();							// generate random number between 0 and 1 
		if (rn2 < pst) {
			con.put(id, nx, ny, nz);		// pridani castice
			std::cout << "YES \n";
			if (id > npart) {				// pouzilo-li se ID o jedna vyssi nez maximalni ID v kontejneru a pridani
				npart = npart + 1;			// castice se uskutecni, pak je nutne navysit toto maximalni ID o 1
			}
			else {							// pouzilo-li se ID z fid, musi ve fid smazat
				fid.erase(fid.begin());
			}
			rn_add++;						// counter of added particles
		}
		else { std::cout << "NO \n"; }
		//std::cout << npart << "\n";  /////////////////////////////////////////////////////////////////////////////////
	}     // BIRTH

	if (((const1) < rn1) && (rn1 <= (const2))) {
		del = uniform_int(1, con.total_particles());	// vybere castici (resp. jeji poradi v containeru) ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &con);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = con.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		std::cout << "DELETE " << del << " " << con.p[ijk][3 * q] << " " << con.p[ijk][3 * q + 1] << " " << con.p[ijk][3 * q + 2] << "\n"; // popis cinosti
														// nyni je potreba tuto castici odstranit z datovych struktur containeru, nestane se container nestabilni???
														// pred vymazanim je treba si zapsat jeji sousedy:

		no = con.total_particles();                     // pocet castic v kontejneru
		pst = try_delete_2(ijk, q, con, theta, alfa, beta, B, fc);		// spocte pravdepodobnost se kterou dojde k operaci DELETE
																		// pst = 0.5;
		std::cout << pst << " ";
		pst = pst*(no / zet);
		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {
			erase(ijk, q, &fid, &con);						// smazani castice - uvolni ID do fid k opetovnemu pouziti
			std::cout << "YES \n";
			rn_del++;										// counter of deleted particles
		}
		else { std::cout << "NO \n"; }
	}     // DEATH

	if (rn1 > (const2)) {
		del = uniform_int(1, con.total_particles());	// vybere castici ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &con);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = con.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = con.p[ijk][3 * q]; y = con.p[ijk][3 * q + 1]; z = con.p[ijk][3 * q + 2]; // urci jeji souradnice
		// vygeneruje novou castici z normalniho rozdeleni:
		nx = normal(x, sigma); ny = normal(y, sigma); nz = normal(z, sigma);   // coordinates of new particle
		//std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
		// these coordinates could be outside the window [0,1]^3 and in this case it is necessary to remap it into window:
		nx = nx - step_int(nx); ny = ny - step_int(ny); nz = nz - step_int(nz);
		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti

		pst = try_move_2(ijk, q, nx, ny, nz, con, theta, alfa, beta, B, fc);		// urci pravdepodobnost se kterou dojde k operaci MOVE
																					// pst = 0.5;
		std::cout << pst << " ";

																					// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {
			erase(ijk, q, &con);					// smazani castice
			con.put(del, nx, ny, nz);				// pridani nove castice (ID zustava zachovano)

			std::cout << "YES \n";
			rn_mov++;								// counter of moved particles
		}
		else { std::cout << "NO \n"; }
	}     // MOVE 
}
*/

// bez vypisovani:
void bdma_step_2(long &npart, std::vector<int> &fid, voro::container &con, double sigma, double theta, double zet, const double &alfa, const double &beta, const double &B, long &rn_add, long &rn_del, long &rn_mov, bool &fc)
{
	// [in,out]		npart	number of added particles in the whole history.
	// [in,out]		fid		vector of available id큦.
	// [in]			con		the container with stored particles.
	// [in]			sigma	parameter of normal distribution.
	// [in]			theta	parameter of V2 function
	// [in]			zet     intensity constant of reference Poisson process
	// [in]			fc		identificator of energy (0 - uses face areas, 1 - uses dihedral angles)

	double rn1 = uniform(0, 1);   // random number between 0 and 1 
	double rn2 = uniform(0, 1);   // random number between 0 and 1 
								  //rn1 = 0.8;
								  //rn2 = 0.0000000000000001;
	double nx, ny, nz, x, y, z;
	double pst;
	double const1 = 1;
	double const2 = 2;
	double const3 = 3;

	int del, id, no;
	int ijk, q; // ijk_add, q_add; unsigned int i,j,k; 

	
	const1 = const1 / const3;
	const2 = const2 / const3;

	if (rn1 <= (const1)) {
		nx = uniform(0, 1); ny = uniform(0, 1); nz = uniform(0, 1); // coordinates of new particle
																	// what about ID ??? - dve moznosti - pamatovat pocet jiz pridanych castic, tj. pamatovat nejvyssi jiz pouzite ID
																	//                                  - vezt zaznam o uvolnenych ID po odstranenych casticich a znovu je pouzivat
		if (fid.size() > 0) {				// najiti vhodneho ID
			id = fid[0];
		}   // zamerem bylo zjistit zdali je ve fid ulozene nejake volne id ci ne???
		else { id = npart + 1; }  // pokud nebylo dostupne ID ve fid, pouzije se maximalni hodnota ID + 1
//		std::cout << "ADD " << id << " " << nx << " " << ny << " " << nz << " \n"; // popis cinosti

		no = con.total_particles();                             // pocet castic v kontejneru
		pst = try_add_2(id, nx, ny, nz, con, theta, alfa, beta, B, fc);		  	// spocte pravdepodobnost se kterou dojde k operaci ADD
																			// pst = 0.5;
																			//		std::cout << pst << " ";
		pst = pst*(zet / (no + 1));
//		std::cout << pst << " ";

		// rn = rnd();							// generate random number between 0 and 1 
		if (rn2 < pst) {
			con.put(id, nx, ny, nz);		// pridani castice
//			std::cout << "YES \n";
			if (id > npart) {				// pouzilo-li se ID o jedna vyssi nez maximalni ID v kontejneru a pridani
				npart = npart + 1;			// castice se uskutecni, pak je nutne navysit toto maximalni ID o 1
			}
			else {							// pouzilo-li se ID z fid, musi ve fid smazat
				fid.erase(fid.begin());
			}
			rn_add++;						// counter of added particles
		}
//		else { std::cout << "NO \n"; }
		//std::cout << npart << "\n";  /////////////////////////////////////////////////////////////////////////////////
	}     // BIRTH

	if (((const1) < rn1) && (rn1 <= (const2))) {
		del = uniform_int(1, con.total_particles());	// vybere castici (resp. jeji poradi v containeru) ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &con);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = con.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
//		std::cout << "DELETE " << del << " " << con.p[ijk][3 * q] << " " << con.p[ijk][3 * q + 1] << " " << con.p[ijk][3 * q + 2] << "\n"; // popis cinosti
														// nyni je potreba tuto castici odstranit z datovych struktur containeru, nestane se container nestabilni???
														// pred vymazanim je treba si zapsat jeji sousedy:

		no = con.total_particles();                     // pocet castic v kontejneru
		pst = try_delete_2(ijk, q, con, theta, alfa, beta, B, fc);		// spocte pravdepodobnost se kterou dojde k operaci DELETE
																	// pst = 0.5;
																	//		std::cout << pst << " ";
		pst = pst*(no / zet);
//		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {
			erase(ijk, q, &fid, &con);						// smazani castice - uvolni ID do fid k opetovnemu pouziti
//			std::cout << "YES \n";
			rn_del++;										// counter of deleted particles
		}
//		else { std::cout << "NO \n"; }
	}     // DEATH

	if (rn1 > (const2)) {
		del = uniform_int(1, con.total_particles());	// vybere castici ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &con);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = con.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = con.p[ijk][3 * q]; y = con.p[ijk][3 * q + 1]; z = con.p[ijk][3 * q + 2]; // urci jeji souradnice
																					 // vygeneruje novou castici z normalniho rozdeleni:
		nx = normal(x, sigma); ny = normal(y, sigma); nz = normal(z, sigma);   // coordinates of new particle
		//std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
																			   // these coordinates could be outside the window [0,1]^3 and in this case it is necessary to remap it into window:
		nx = nx - step_int(nx); ny = ny - step_int(ny); nz = nz - step_int(nz);
//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti

		pst = try_move_2(ijk, q, nx, ny, nz, con, theta, alfa, beta, B, fc);		// urci pravdepodobnost se kterou dojde k operaci MOVE
		// pst = 0.5;
//		std::cout << pst << " ";

																			
		if (rn2 <= pst) {
			erase(ijk, q, &con);					// smazani castice
			con.put(del, nx, ny, nz);				// pridani nove castice (ID zustava zachovano)

//			std::cout << "YES \n";
			rn_mov++;								// counter of moved particles
		}
//		else { std::cout << "NO \n"; }
	}     // MOVE 
}


// bez vypisovani :
void bdma_step_3(long &npart, std::vector<int> &fid, voro::container &con, double sigma, double theta1, double theta2, double zet, const double &alfa, const double &beta, const double &B, long &rn_add, long &rn_del, long &rn_mov, bool &fc)
{
	// [in,out]		npart	number of added particles in the whole history.
	// [in,out]		fid		vector of available id큦.
	// [in]			con		the container with stored particles.
	// [in]			sigma	parameter of normal distribution.
	// [in]			theta	parameter of V2 function
	// [in]			zet     intensity constant of reference Poisson process
	// [in]			fc		identificator of energy (0 - uses face areas, 1 - uses dihedral angles)

	double rn1 = uniform(0, 1);   // random number between 0 and 1 
	double rn2 = uniform(0, 1);   // random number between 0 and 1 
								  //rn1 = 0.8;
								  //rn2 = 0.0000000000000001;
	double nx, ny, nz, x, y, z;
	double pst;
	double const1 = 1;
	double const2 = 2;
	double const3 = 3;

	int del, id, no;
	int ijk, q; // ijk_add, q_add; unsigned int i,j,k; 


	const1 = const1 / const3;
	const2 = const2 / const3;

	if (rn1 <= (const1)) {
		nx = uniform(0, 1); ny = uniform(0, 1); nz = uniform(0, 1); // coordinates of new particle
																	// what about ID ??? - dve moznosti - pamatovat pocet jiz pridanych castic, tj. pamatovat nejvyssi jiz pouzite ID
																	//                                  - vezt zaznam o uvolnenych ID po odstranenych casticich a znovu je pouzivat
		if (fid.size() > 0) {				// najiti vhodneho ID
			id = fid[0];
		}   // zamerem bylo zjistit zdali je ve fid ulozene nejake volne id ci ne???
		else { id = npart + 1; }  // pokud nebylo dostupne ID ve fid, pouzije se maximalni hodnota ID + 1
//		std::cout << "ADD " << id << " " << nx << " " << ny << " " << nz << " \n"; // popis cinosti

		no = con.total_particles();                             // pocet castic v kontejneru

		// kombinace perove a terciarni energie
		pst = try_add_2(id, nx, ny, nz, con, theta2, alfa, beta, B, fc);		  	// spocte pravdepodobnost se kterou dojde k operaci ADD
		pst = pst + try_add(id, nx, ny, nz, con, theta1, alfa, beta, B);

//		std::cout << pst << " ";
		pst = pst*(zet / (no + 1));
//		std::cout << pst << " ";

		// rn = rnd();							// generate random number between 0 and 1 
		if (rn2 < pst) {
			con.put(id, nx, ny, nz);		// pridani castice
//			std::cout << "YES \n";
			if (id > npart) {				// pouzilo-li se ID o jedna vyssi nez maximalni ID v kontejneru a pridani
				npart = npart + 1;			// castice se uskutecni, pak je nutne navysit toto maximalni ID o 1
			}
			else {							// pouzilo-li se ID z fid, musi ve fid smazat
				fid.erase(fid.begin());
			}
			rn_add++;						// counter of added particles
		}
//		else { std::cout << "NO \n"; }
		//std::cout << npart << "\n";  /////////////////////////////////////////////////////////////////////////////////
	}     // BIRTH

	if (((const1) < rn1) && (rn1 <= (const2))) {
		del = uniform_int(1, con.total_particles());	// vybere castici (resp. jeji poradi v containeru) ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &con);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = con.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
//		std::cout << "DELETE " << del << " " << con.p[ijk][3 * q] << " " << con.p[ijk][3 * q + 1] << " " << con.p[ijk][3 * q + 2] << "\n"; // popis cinosti
														// nyni je potreba tuto castici odstranit z datovych struktur containeru, nestane se container nestabilni???
														// pred vymazanim je treba si zapsat jeji sousedy:

		no = con.total_particles();                     // pocet castic v kontejneru

		// kombinace parove a terciarni energie
		pst = try_delete_2(ijk, q, con, theta2, alfa, beta, B, fc);		// spocte pravdepodobnost se kterou dojde k operaci DELETE
		pst = pst + try_delete(ijk, q, con, theta1, alfa, beta, B);

//		std::cout << pst << " ";
		pst = pst*(no / zet);
//		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {
			erase(ijk, q, &fid, &con);						// smazani castice - uvolni ID do fid k opetovnemu pouziti
//			std::cout << "YES \n";
			rn_del++;										// counter of deleted particles
		}
//		else { std::cout << "NO \n"; }
	}     // DEATH

	if (rn1 > (const2)) {
		del = uniform_int(1, con.total_particles());	// vybere castici ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &con);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = con.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = con.p[ijk][3 * q]; y = con.p[ijk][3 * q + 1]; z = con.p[ijk][3 * q + 2]; // urci jeji souradnice
																					 // vygeneruje novou castici z normalniho rozdeleni:
		nx = normal(x, sigma); ny = normal(y, sigma); nz = normal(z, sigma);   // coordinates of new particle
																			   //std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
																			   // these coordinates could be outside the window [0,1]^3 and in this case it is necessary to remap it into window:
		nx = nx - step_int(nx); ny = ny - step_int(ny); nz = nz - step_int(nz);
//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti

		pst = try_move_2(ijk, q, nx, ny, nz, con, theta2, alfa, beta, B, fc);		// urci pravdepodobnost se kterou dojde k operaci MOVE
		pst = pst + try_move(ijk, q, nx, ny, nz, con, theta1, alfa, beta, B);

//		std::cout << pst << " ";

		if (rn2 <= pst) {
			erase(ijk, q, &con);					// smazani castice
			con.put(del, nx, ny, nz);				// pridani nove castice (ID zustava zachovano)

//			std::cout << "YES \n";
			rn_mov++;								// counter of moved particles
		}
//		else { std::cout << "NO \n"; }
	}     // MOVE 
}


// fc V1 computes mean number of faces for the primary and the secondary particles
double V1(voro::container &con, int pr, std::vector<int> &rsr)
{
	//	[in]	con		the container with stored particles.
	//	[in]	pr		number of faces of primary particle, if 0 then there is no primary particle
	//	[in]	sr		list of secondary particles

	voronoicell_neighbor c;
	double mean = pr;
	int counter = 0;
	unsigned int i;

	if (pr > 0) { counter++; }										// increase counter in case there was a primary particle

	for (i = 0; i < rsr.size() / 2; i++) {                          // loop over secondary particles
		// ijk_i = sr[2 * i]; q_i = sr[2 * i + 1];
		con.compute_cell(c, rsr[2 * i], rsr[2 * i + 1]);

		mean = mean + c.number_of_faces();
		counter++;
	}


	return (mean/counter);
}


// fc V1 computes energy of a given cell with respect to number of neighbours (capture only the mean number not variance)
// cell with similar number of neighbours to the mean number of neighbours computed from data is preferred
double V1_nf(voronoicell_neighbor &rc, int &n0)
{
	// [in]		c				the cell for which is function V1_nf computed.
	// [in]		n0				preferable number of neighbours

	int n;
	int m = 2; // exponent

	n = rc.number_of_faces();

	return ( abs(n - n0) / n0) ^ m;
}


// fce V2 pro pocitani podilu objemu 
// pouzit jako argumenty primo bunky a nebo pozice danych castic (ijk,q) ?; pokud bunky, budu pak potrebovat informaci
//   o sousedech ? (kdyby ano musim pouzit voronoicell_neighbor misto voronoicell)
double V2(voronoicell_neighbor &rc1, voronoicell_neighbor &rc2, double &rx, double &ry, double &rz, double &rxx, double &ryy, double &rzz)
{
	// [in]		c1,c2				the cells for which is function V2 computed.
	// [in]		x,y,z,xx,yy,zz		the coordinates of cell generators

	double x1, y1, z1, x2, y2, z2;

	double a = rc1.volume();  // vypocet objemu
	double b = rc2.volume();
	//std::cout << "Volumes: " << a << " " << b << "\n";

	rc1.centroid(x1, y1, z1);	// vypocet tezist obou bunek (v relativnich souradnicich vzhledem ke generatoru!!!)
	rc2.centroid(x2, y2, z2);
	
	x1 = rx + x1; y1 = ry + y1; z1 = rz + z1;		// skutecne souradnice tezist
	x2 = rxx + x2; y2 = ryy + y2; z2 = rzz + z2;
	//std::cout << "Centroids: (1) " << x1 << " " << y1 << " " << z1 << " (2) " << x2 << " " << y2 << " " << z2 << "\n";

	if (barycentrum(x1, y1, z1, x2, y2, z2, a, b)) {	// je-li teziste sjednoceni bunek v okne spocti hodnotu, jinak vrat 0

		min_max(a, b);		// a >= b

		// !!!! POZOR konstanta pro porovnani vychazi ze vzorce pro objem koule (4/3 Pi r^3) coz zavisi na hodnote r = alfa
		// a se neporovnava prtz a>b
		// pro alfa=0.05 je minimalni objem 0.00006544984695
		//if (b < 0.00006544984695) { std::cout << "PODEZRELE MALY OBJEM BUNKY! - " << b << " ; " << a << " \n"; /*std::cin.ignore();*/ };
		// konstanta spoctena pro alfa = 0.05, a tedy min polomer bunky alfa/2 = 0.025

		return sqrt(a / b - 1);
		// existuji fce std::max a std::min
	}
	else { return 0; }
	
}



// fc V3f returns V3 value (wrt face areas) of given three cells with common edge
// VSTUPEM MUSI BYT TRI SOUSEDICI BUNKY SDILEJICI HRANU !!!!!!!!!!!!!!!! - overit
// - nestaci tri sousedici bunky, nebot i kdyz sousedi, tak jeste nemusi urcovat hranu ??? - hypoteza
// - klicove je tedy overit, zda tri bunky urcuji hranu, pak lze predpokladat, ze i navzajem sousedi
// VSTUPEM JSOU TRI NAVZAJEM SOUSEDICI BUNKY, JE JESTE TREBA OVERIT, ZDA URCUJI HRANU
double V3f(int &id1, voronoicell_neighbor &rc1, int &id2, voronoicell_neighbor &rc2, int &id3, voronoicell_neighbor &rc3,
				double &rx, double &ry, double &rz, double &rxx, double &ryy, double &rzz, double &rxxx, double &ryyy, double &rzzz)
{
	// [in]		id1, id2, id3		ID numbers of three cells sharing the edge, ordered by value, i.e. id1 < id2 < id3
	// [in]		c1, c2, c3			the three cells sharing an edge
	// [in]		x,y,z,xx, ...		the true coordinates of generators of these cells

//	int i;
	int i1, i2, i3;
	double a1,a2,a3; 
	std::vector<int> neigh;
	std::vector<double> areas;

	double x1, y1, z1, x2, y2, z2;

	
	// hrana je urcena tremi bunkami. Postup:
// idea) pro prvni bunku najdi steny prislusne druhe a treti bunce. Pro druhou najdi stenu prislusnou treti bunce. Tj. mame vsechny tri steny
//	sdilejici zadanou hranu. Spocteme tedy jejich obsahy a jsme hotovi.

		rc1.neighbors(neigh);
		i1 = find_index(neigh, id2);
		if (i1 == -1) { return 0; }
		//if (i1 == -1) { std::cout << "V3f (1) - not neighbouring particles:" << id2 << " \n"; return 0; }
		//for (i = 0; i < neigh.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
		//	std::cout << neigh[i] << " ";
		//} std::cout << " \n"; return 0; }
		i2 = find_index(neigh, id3);
		if (i2 == -1) { return 0; }
		//if (i2 == -1) { std::cout << "V3f (2) - not neighbouring particles:" << id3 << " \n"; return 0; }

// 1) over zda trojice bunek urcuje stenu
		if (common_edge(rc1, i1, i2) == -1) { /*std::cout << "V3f - no common edge \n";*/ return 0; }
		// overeni podminky, zda vubec dana trojice bunek urcuje hranu (tj. zda-li prispiva nejakou energii)
		// ??? - je lepsi prvne overit polohu teziste, nebo zda urcuji hranu (u obou poruseni vede k navratu nulove energie) - co je castejsi?
		// odpoved: (odkomentovanim zjistis cetnost) mnohem castejsi nez poloha teziste mimo okno (nastava jen u hranic, kdezto problem urceni hrany je celoplosny)

// 2) korekce vuci periodicite muze byt jako v pripade V2, tj uvazovat teziste trojice bunek - pak ale musime uvazovat terciarni strukturu i ve fci V3a;
//		druhou moznosti je testovat, zda stred dotycne hrany je v okne ci nikoliv - umoznuje zjednodusit vypocet, nebot pri pouziti V3a se obejdeme bez terciarni struktury
// 2a) over zda teziste trojice je uvnitr okna
		double a = rc1.volume();  // vypocet objemu
		double b = rc2.volume();
		double c = rc3.volume();

		rc1.centroid(x1, y1, z1);	// vypocet tezist prvnich dvou bunek (v relativnich souradnicich vzhledem ke generatoru!!!)
		rc2.centroid(x2, y2, z2);

		x1 = rx + x1; y1 = ry + y1; z1 = rz + z1;		// skutecne souradnice tezist
		x2 = rxx + x2; y2 = ryy + y2; z2 = rzz + z2;

		bar_coor(x1, y1, z1, x2, y2, z2, a, b); // spocte spolecne teziste prvnich dvou bunek do vektoru (x1,y1,z1)
		rc3.centroid(x2, y2, z2);
		x2 = rxxx + x2; y2 = ryyy + y2; z2 = rzzz + z2;

		if (barycentrum(x1, y1, z1, x2, y2, z2, a + b, c) == false) {
			return 0;
		}	// je-li teziste sjednoceni bunek v okne spocti hodnotu, jinak vrat 0

// 2b) over zda stred hrany je uvnitr okna
// --- zde neni nutne potreba


		rc2.neighbors(neigh);
		i3 = find_index(neigh, id3);
		if (i3 == -1) { std::cout << "ERROR: V3f (3) - not neighbouring particles \n"; return 0; }

// 3) urci energii
		rc1.face_areas(areas);
		a1 = areas[i1];
		a2 = areas[i2];
		rc2.face_areas(areas);
		a3 = areas[i3];

		min_max(a1, a2, a3);

		return sqrt(a3 / a1 - 1);
	//}
	//else { return 0; }
}

// fc V3a returns V3 value (wrt dihedral angles) of given three cells with common edge
// VSTUPEM MUSI BYT TRI SOUSEDICI BUNKY SDILEJICI HRANU !!!!!!!!!!!!!!!! - overit
// - nestaci tri sousedici bunky, nebot i kdyz sousedi, tak jeste nemusi urcovat hranu ??? - hypoteza
// - klicove je tedy overit, zda tri bunky urcuji hranu, pak lze predpokladat, ze i navzajem sousedi
// VSTUPEM JSOU TRI NAVZAJEM SOUSEDICI BUNKY, JE JESTE TREBA OVERIT, ZDA URCUJI HRANU
double V3a(int &id1, voronoicell_neighbor &rc1, int &id2, voronoicell_neighbor &rc2, int &id3, voronoicell_neighbor &rc3,
				double &rx, double &ry, double &rz, double &rxx, double &ryy, double &rzz, double &rxxx, double &ryyy, double &rzzz)
{
	// [in]		id1, id2, id3		ID numbers of three cells sharing the edge, ordered by value, i.e. id1 < id2 < id3
	// [in]		c1, c2, c3			the first, the second and the third cell sharing an edge
	// [in]		x,y,z,xx, ...		the true coordinates of generators of these cells

	int i1, i2, i, j, k, p;
	double a1, a2, a3;
	std::vector<int> neigh;
	std::vector<double> angles;

	double x1, y1, z1;

// 1a) v prvni bunce urcim steny prislusne druhe a treti bunce (sousedi-li bunky). Tyto steny sdili jednu hranu, tu najdu (tj. najdu jeji dva vrcholy a urcim 
//	poradi teto hrany ve strukture ed prvni bunky).  (pokud takova hrana vubec existuje)
		rc1.neighbors(neigh);
		i1 = find_index(neigh, id2);
		if (i1 == -1) { return 0; }
		//if (i1 == -1) { std::cout << "V3a (1) - not neighbouring particles \n"; return 0; }
		i2 = find_index(neigh, id3);
		if (i2 == -1) { return 0; }
		//if (i2 == -1) { std::cout << "V3a (1) - not neighbouring particles \n"; return 0; }

		// i1, i2 jsou nyni cisla sten

		// pro dve steny chci najit spolecnou hranu
		p = common_edge(rc1, i1, i2);		
		if (p == -1) { return 0; }
		//if (p == -1) { std::cout << "V3a (1) - no common edge \n"; return 0; }

		// pozn.: pokud se v teto prvni casti nevypise "not neigh particles" nebo "no common edge", pak uz k tomu nedojde ani v druhe ci treti casti

// 2) korekce vuci periodicite muze byt jako v pripade V2, tj uvazovat teziste trojice bunek - pak ale musime uvazovat terciarni strukturu i ve fci V3a;
//		druhou moznosti je testovat, zda stred dotycne hrany je v okne ci nikoliv - umoznuje zjednodusit vypocet, nebot pri pouziti V3a se obejdeme bez terciarni struktury
// 2a) over zda teziste trojice je uvnitr okna
/*		double a = rc1.volume();  // vypocet objemu
		double b = rc2.volume();
		double c = rc3.volume();

		rc1.centroid(x1, y1, z1);	// vypocet tezist prvnich dvou bunek (v relativnich souradnicich vzhledem ke generatoru!!!)
		rc2.centroid(x2, y2, z2);

		x1 = rx + x1; y1 = ry + y1; z1 = rz + z1;		// skutecne souradnice tezist
		x2 = rxx + x2; y2 = ryy + y2; z2 = rzz + z2;

		bar_coor(x1, y1, z1, x2, y2, z2, a, b); // spocte spolecne teziste prvnich dvou bunek do vektoru (x1,y1,z1)
		rc3.centroid(x2, y2, z2);
		x2 = rxxx + x2; y2 = ryyy + y2; z2 = rzzz + z2;

		if (barycentrum(x1, y1, z1, x2, y2, z2, a + b, c) == false) {
			return 0;
		}	// je-li teziste sjednoceni bunek v okne spocti hodnotu, jinak vrat 0
*/ 
// 2b) over zda stred hrany je uvnitr okna
		// p - zkoumana hrana - v bunce c1
		i = 0; j = 0; k = 0;
		//number_of_vertices = 2 + rc1.number_of_edges() - rc1.number_of_faces(); // porovnat s rc1.p
		//std::cout << rc1.p << "  vs  " << number_of_vertices << "\n";

//		int v1, v2;
		while (i < p+1) {
			if (k < rc1.nu[j]) {					// loop pres hrany: c.p = snad pocet vrcholu, c.nu - vektor radu vrcholu
				if (j < rc1.ed[j][k]) { i++; }
				k++;
			}
			else { k = 0; j++; }	
		}	// vysledkem je, ze j a ed[j][k-1] jsou vrcholy hledane hrany

		// x1,y1,z1 - souradnice stredu strany
		x1 = (rx + rc1.pts[3 * j] * 0.5) + (rx + rc1.pts[3 * rc1.ed[j][k - 1]] * 0.5); x1 = x1 / 2;
		y1 = (ry + rc1.pts[3 * j + 1] * 0.5) + (ry + rc1.pts[3 * rc1.ed[j][k - 1] + 1] * 0.5); y1 = y1 / 2;
		z1 = (rz + rc1.pts[3 * j + 2] * 0.5) + (rz + rc1.pts[3 * rc1.ed[j][k - 1] + 2] * 0.5); z1 = z1 / 2;

		if (x1 > 0 && x1 < 1 && y1 > 0 && y1 < 1 && z1 > 0 && z1 < 1) {}
		else { return 0; }

// 3) Znam-li poradi hrany, mohu spocitat jeji dihedralni uhel prislusny prvni bunce.
		dihedral_angles(angles, rc1);
		a1 = angles[p];					// dihedralni uhel v radianech mezi stenami prvni bunky

// 4) Abych urcil i zbyvajici dva dihedralni uhly, musim kroky 1, 2 a 3 zopakovat i pro zbyle dve bunky.
		rc2.neighbors(neigh);
		i1 = find_index(neigh, id1);
		if (i1 == -1) { std::cout << "ERROR: V3a (2) - not neighbouring particles \n"; return 0; }
		i2 = find_index(neigh, id3);
		if (i2 == -1) { std::cout << "ERROR: V3a (2) - not neighbouring particles \n"; return 0; }

		p = common_edge(rc2, i1, i2);
		if (p == -1) { std::cout << "ERROR: V3a (2) - no common edge \n"; return 0; }

		dihedral_angles(angles, rc2);
		a2 = angles[p];					// dihedralni uhel v radianech mezi stenami druhe bunky

		rc3.neighbors(neigh);
		i1 = find_index(neigh, id1);
		if (i1 == -1) { std::cout << "ERROR: V3a (3) - not neighbouring particles \n"; return 0; }
		i2 = find_index(neigh, id2);
		if (i2 == -1) { std::cout << "ERROR: V3a (3) - not neighbouring particles \n"; return 0; }

		p = common_edge(rc3, i1, i2);
		if (p == -1) { std::cout << "ERROR: V3a (3) - no common edge \n"; return 0; }

		dihedral_angles(angles, rc3);
		a3 = angles[p];					// dihedralni uhel v radianech mezi stenami treti bunky


		min_max(a1, a2, a3);

		return sqrt(a3 / a1 - 1);
	
}

// observation: the vertex is determined by four cells
double V4(int &id1, voronoicell_neighbor &rc1, int &id2, voronoicell_neighbor &rc2, int &id3, int &id4)
{
	// [in]		id1, id2, id3, id4		ID numbers of four cells sharing the vertex, ordered by value, i.e. id1 < id2 < id3 < id4
	// [in]		c1, c2					the first and the second cell

	int i1, i2, i3, p1, p2, p3, v1, v2;
	int k, l, m;
	double l1, l2;
	std::vector<int> neigh, ord, vert;
	std::vector<double> lengths, le;
	lengths.resize(4);

	// ctyri bunky urcuji vrchol. Postup:
// 1) pro prvni bunku najdu pozice zbylych tri sousednich bunek, tim najdu tri steny. 
	
	rc1.neighbors(neigh);
	i1 = find_index(neigh, id2);	// v prvni bunce najdu tri steny prislusne trem sousednim bunkam
	i2 = find_index(neigh, id3);	// tyto steny budou mit jeden spolecny vrchol
	i3 = find_index(neigh, id4);

	rc1.face_orders(ord);		// poznacim si vzdy prvni vrchol techto sten ve vektoru face_vertices
	p1 = 0;
	for (k = 0; k < i1; k++) {
		p1 = p1 + ord[k];
	}
	p1 = p1 + i1 + 1;			// p1 - prvni vrchol steny i1
	p2 = 0;
	for (k = 0; k < i2; k++) {
		p2 = p2 + ord[k];
	}
	p2 = p2 + i2 + 1;			// p2 - prvni vrchol steny i2
	p3 = 0;
	for (k = 0; k < i3; k++) {
		p3 = p3 + ord[k];
	}
	p3 = p3 + i3 + 1;			// p3 - prvni vrchol steny i3

// 2) Tyto tri steny maji spolecny prave jeden vrchol. Najdu ho.
	rc1.face_vertices(vert);
	for (k = 0; k < ord[i1]; k++) {				//najdi spolecny vrchol techto sten
		for (l = 0; l < ord[i2]; l++) {
			for (m = 0; m < ord[i3]; m++) {
				if (vert[p1 + k] == vert[p2 + l] && vert[p2 + l] == vert[p3 + m]) { v1 = vert[p1 + k]; break; }
			}
		}
	}

// 3) Pro tento vrchol urcim jeho sousedni vrcholy, dostanu tak tri hrany teto bunky. 
// 4) Urcim poradi techto hran ve strukture ed a spoctu delky techto tri hran.
	edge_lengths(le, rc1);
	for (k = 0; k < rc1.nu[v1]; k++) {
		v2 = rc1.ed[v1][k];
		min_max(v1, v2);			// v1 >= v2

		m = 0;										// najdi poradi teto hrany
		for (l = 0; l < v2; l++) {
			m = m + rc1.nu[l];
		}
		l = 0;
		m++;
		while (rc1.ed[v2][l] != v1) {
			m++; l++; if (l > rc1.nu[v2]) { break; std::cout << " ERROR: Edge not found! \n"; }
		}
		lengths[k] = le[m];
	}

// zbyva posledni ctvrta hrana, ta uz je ale soucasti jine bunky (kterekoliv ze tri sousednich):
// 5) pro druhou bunku najdu pozice ostatnich tri bunek, a tedy prislusnych sten, i1 bude pozice prvni bunky. Najdu pozice prvnich vrcholu 
//	techto sten ve vektoru face_vertices. Urcim jejich spolecny bod (tentyz bod co v bodu 2, jen ted muze mit jine poradove cislo).
	rc2.neighbors(neigh);
	i1 = find_index(neigh, id1);
	i2 = find_index(neigh, id3);	// tyto steny budou mit opet jeden spolecny vrchol
	i3 = find_index(neigh, id4);

	rc2.face_orders(ord);
	p1 = 0;
	for (k = 0; k < i1; k++) {
		p1 = p1 + ord[k];
	}
	p1 = p1 + i1 + 1;			// p1 - prvni vrchol steny i1
	p2 = 0;
	for (k = 0; k < i2; k++) {
		p2 = p2 + ord[k];
	}
	p2 = p2 + i2 + 1;			// p2 - prvni vrchol steny i2
	p3 = 0;
	for (k = 0; k < i3; k++) {
		p3 = p3 + ord[k];
	}
	p3 = p3 + i3 + 1;			// p3 - prvni vrchol steny i3

	rc2.face_vertices(vert);
	for (k = 0; k < ord[i1]; k++) {				//najdi spolecny vrchol techto sten
		for (l = 0; l < ord[i2]; l++) {
			for (m = 0; m < ord[i3]; m++) {
				if (vert[p1 + k] == vert[p2 + l] && vert[p2 + l] == vert[p3 + m]) { v1 = vert[p1 + k]; break; }
			}
		}
	}

// 6) Ze tri sousedu tohoto vrcholu budou dva lezet v prvni stene a treti nikoliv. Ten najdu.
	for (k = 0; k < rc2.nu[v1]; k++) {			// najdi sousedni vrcchol, ktery nelezi v prvni stene
		v2 = rc2.ed[v1][k];
		m = 0;
		for (l = 0; l < ord[i1]; l++) {		
			if (v2 == vert[p1 + l]) { m = 1; }  // lezi v prvni stene, pokracuj
		}
		if (m == 0) { break; }					// nasli jsme vrchol, ktery neni v prvni stene, ukonci cyklus
	}

// 7) Pro dva nalezene vrcholy, tj. pro posledni hranu, urci jeji poradi ve strukture ed a spocti jeji delku.
	min_max(v1, v2);			// v1 >= v2

	m = 0;										// najdi poradi teto hrany
	for (l = 0; l < v2; l++) {
		m = m + rc2.nu[l];
	}
	l = 0;
	m++;
	while (rc2.ed[v2][l] != v1) {
		m++; l++; if (l > rc2.nu[v2]) { break; std::cout << " ERROR: Edge not found! \n"; }
	}
	edge_lengths(le, rc2);
	lengths[4] = le[m];

// HOTOVO, mam delky vsech ctyr hran vybihajicich z urceneho vrcholu
	l1 = 1; l2 = 0;
	for (k = 0; k < lengths.size(); k++) {
		if (lengths[k] < v1) { l1 = lengths[k]; }
		if (lengths[k] > v2) { l2 = lengths[k]; }
	}

	return sqrt(l2 / l1 - 1);
}



// fce BDMA_step - telo algoritmu, provadi operace ADD, DELETE, MOVE, CHANGE of RADIUS a vola prislusne funkce TRY pro spocteni prislusnych psti; radical version
void bdma_step(long &npart, std::vector<int> &fid, voro::container_poly &conp, double sigma, double theta, double zet, const double &alfa, const double &beta, const double &B, const double &iota, long &rn_add, long &rn_del, long &rn_mov, double &arad)
{
	// [in,out]		npart						number of added particles in the whole history.
	// [in,out]		fid							vector of available id큦.
	// [in]			con							the container with stored particles.
	// [in]			sigma						parameter of normal distribution.
	// [in]			theta						parameter of V2 function
	// [in]			zet							intensity constant of reference Poisson process
	// [in]			alfa, beta, B				hardcore parameters
	// [out]		n_add, n_del, n_mov, n_chr	counters of added/deleted/moved/radius changed particles
	// [in]			arad						specifies the distribution of radii - zatim nadbytecne

	double rn1 = uniform(0, 1);   // random number between 0 and 1 
	double rn2 = uniform(0, 1);   // random number between 0 and 1 
	//rn1 = 0.2;
	//rn2 = 0.0000000000000001;
	double nx, ny, nz, x, y, z, r, nr;
	double pst;
	double const1 = 1;
	double const2 = 2;
	double const3 = 3;

	std::cout << "jsem tu";
/*
	int del, id, no;
	int ijk, q; // ijk_add, q_add; unsigned int i,j,k; 

	//std::cout << "pst: " << rn1 << " " << rn2 << "\n";
	const1 = const1 / const3;
	const2 = const2 / const3;

	// std::cout << "fid: ";	///////////////////////////////////////////////////////////////////////////////////////////
	// for (i = 0; i < fid.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	// std::cout << fid[i] << " ";
	// } std::cout << " \n";

	if (rn1 <= (const1)) {
		nx = uniform(0, 1); ny = uniform(0, 1); nz = uniform(0, 1); // coordinates of new particle
		// the radius of the new particle can be generated from the prespecified distribution:
		//r = 0.005 + 0.025*uniform(0, 1); 
		r = 0.05*uniform(0, 1);
		//r = gamma(3, 0.5);
		//r = triangle(0.005, 0.03, 0.02);
		// or using the average radius in container:
		//r = ave_rad(conp);
																
		if (fid.size() > 0) {				// najiti vhodneho ID
			id = fid[0];
		}   // zamerem bylo zjistit zdali je ve fid ulozene nejake volne id ci ne???
		else { id = npart + 1; }  // pokud nebylo dostupne ID ve fid, pouzije se maximalni hodnota ID + 1
//		std::cout << "ADD " << id << " " << nx << " " << ny << " " << nz << " " << r << " \n"; // popis cinosti

		no = conp.total_particles();                             // pocet castic v kontejneru
		pst = try_add(id, nx, ny, nz, r, conp, theta, alfa, beta, B, iota);		  	// spocte pravdepodobnost se kterou dojde k operaci ADD
																			// pst = 0.5;
//		std::cout << pst << " ";
		pst = pst*(zet / (no + 1));
//		std::cout << pst << " ";

		// rn = rnd();							// generate random number between 0 and 1 
		if (rn2 < pst) {
			conp.put(id, nx, ny, nz, r);		// pridani castice
//			std::cout << "YES \n";
			if (id > npart) {				// pouzilo-li se ID o jedna vyssi nez maximalni ID v kontejneru a pridani
				npart = npart + 1;			// castice se uskutecni, pak je nutne navysit toto maximalni ID o 1
			}
			else {							// pouzilo-li se ID z fid, musi ve fid smazat
				fid.erase(fid.begin());
			}
			rn_add++;						// counter of added particles
		}
//		else { std::cout << "NO \n"; }
		//		std::cout << npart << "\n";  /////////////////////////////////////////////////////////////////////////////////
	}     // BIRTH

	if (((const1) < rn1) && (rn1 <= (const2))) {
		del = uniform_int(1, conp.total_particles());	// vybere castici (resp. jeji poradi v containeru) ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
//		std::cout << "DELETE " << del << " " << conp.p[ijk][3 * q] << " " << conp.p[ijk][3 * q + 1] << " " << conp.p[ijk][3 * q + 2] << " " << conp.p[ijk][4 * q + 3] << "\n"; // popis cinosti
																																		   // nyni je potreba tuto castici odstranit z datovych struktur containeru, nestane se container nestabilni???
																																		   // pred vymazanim je treba si zapsat jeji sousedy:

		no = conp.total_particles();                     // pocet castic v kontejneru
		pst = try_delete(ijk, q, conp, theta, alfa, beta, B, iota);		// spocte pravdepodobnost se kterou dojde k operaci DELETE
																	// pst = 0.5;
//		std::cout << pst << " ";
		pst = pst*(no / zet);
//		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {
			erase(ijk, q, &fid, &conp);						// smazani castice - uvolni ID do fid k opetovnemu pouziti
//			std::cout << "YES \n";
			rn_del++;										// counter of deleted particles
		}
//		else { std::cout << "NO \n"; }
	}     // DEATH

	if (((const2) < rn1) ) {
		del = uniform_int(1, conp.total_particles());	// vybere castici ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; // urci jeji souradnice
		r = conp.p[ijk][4 * q + 3];	// a jeji radius
		// vygeneruje novou castici z normalniho rozdeleni:
		nx = normal(x, sigma); ny = normal(y, sigma); nz = normal(z, sigma);   // coordinates of new particle
		//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
		// these coordinates could be outside the window [0,1]^3 and in this case it is necessary to remap it into window:
		nx = nx - step_int(nx); ny = ny - step_int(ny); nz = nz - step_int(nz);
		nr = 0.05*uniform(0, 1);
//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti

		pst = try_MOVE(ijk, q, nx, ny, nz, nr, conp, theta, alfa, beta, B, iota);		// urci pravdepodobnost se kterou dojde k operaci MOVE
																			// pst = 0.5;
//		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {
			erase(ijk, q, &conp);					// smazani castice
			conp.put(del, nx, ny, nz, nr);				// pridani nove castice (ID zustava zachovano)

//			std::cout << "YES \n";
			rn_mov++;								// counter of moved particles
		}
//		else { std::cout << "NO \n"; }
	}     // MOVE 
*/

/*	if (((const3) < rn1) ) {
		del = uniform_int(1, conp.total_particles());	// vybere castici ke zmene (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		//x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; // urci jeji souradnice
		r = conp.p[ijk][4 * q + 3];	// urci jeji polomer

		// vygeneruje novy polomer:
		// 1) z useknuteho normalniho rozdeleni:
		//nr = normal(r, 0.0001);	// radius of new particle
		// 2) z exponencialniho rozdeleni
		//nr = exponential(50);			// double lambda = 50;
		// 3) z trohuhelnikoveho rozdeleni
		nr = triangle(0.005, 0.03, 0.02);	// a = 0.005, b = 0.03, c = 0.02
		// pozn.: exponencialni rozdeleni dava nezaporne hodnoty, trojuhelnikove navic i omezene; zato normalni rozdeleni je na cele realne ose, pro 
		//	zjednoduseni ale uvazujme pouze jeho kladnou cast
		//nr = 0.005 + 0.025*uniform(0, 1); 
		//nr = 0.05*uniform(0, 1);
		//nr = gamma(3, 0.5);
		// or using the average radius in container:
		//nr = ave_rad(conp);

//		std::cout << "CHANGE of RADIUS " << del << " " << r << " --> " << nr << "\n"; // popis cinosti
		if (nr < 0) {  } // not allowing negative radius
		else {
			pst = try_change_rad(ijk, q, nr, conp, theta, alfa, beta, B, iota);		// urci pravdepodobnost se kterou dojde k operaci MOVE
																				// pst = 0.5;
//			std::cout << pst << " ";

			// rn = rnd();  // generate random number between 0 and 1 
			if (rn2 <= pst) {
				conp.p[ijk][4 * q + 3] = nr;

//				std::cout << "YES \n";
				rn_chr++;								// counter of particles with changed radii
			}
//			else { std::cout << "NO \n"; }
		}

	}	  // CHANGE of RADIUS*/
}

/*	
void bdma_step_2(long &npart, std::vector<int> &fid, voro::container_poly &conp, double sigma, double theta, double zet, const double &alfa, const double &beta, const double &B, long &rn_add, long &rn_del, long &rn_mov, long &rn_chr, bool &fc)
{
	// [in,out]		npart	number of added particles in the whole history.
	// [in,out]		fid		vector of available id큦.
	// [in]			con		the container with stored particles.
	// [in]			sigma	parameter of normal distribution.
	// [in]			theta	parameter of V2 function
	// [in]			zet     intensity constant of reference Poisson process
	// [in]			fc		identificator of energy (0 - uses face areas, 1 - uses dihedral angles)

	double rn1 = uniform(0, 1);   // random number between 0 and 1 
	double rn2 = uniform(0, 1);   // random number between 0 and 1 
								  //rn1 = 0.8; 
								  //rn2 = 0.0000000000000001;
	double nx, ny, nz, x, y, z, r, nr;
	double pst;
	double const1 = 1;
	double const2 = 2;
	double const3 = 3;
	double const4 = 4;

	int del, id, no;
	int ijk, q; // ijk_add, q_add; unsigned int i,j,k; 


				//std::cout << "pst: " << rn1 << " " << rn2 << "\n";
	const1 = const1 / const4;
	const2 = const2 / const4;
	const3 = const3 / const4;

	if (rn1 <= (const1)) {
		nx = uniform(0, 1); ny = uniform(0, 1); nz = uniform(0, 1); // coordinates of new particle
		r = 0.005 + 0.025*uniform(0, 1);
		r = triangle(0.005, 0.03, 0.02);
		
		if (fid.size() > 0) {				// najiti vhodneho ID
			id = fid[0];
		}   // zamerem bylo zjistit zdali je ve fid ulozene nejake volne id ci ne???
		else { id = npart + 1; }  // pokud nebylo dostupne ID ve fid, pouzije se maximalni hodnota ID + 1
//		std::cout << "ADD " << id << " " << nx << " " << ny << " " << nz << " " << r << " \n"; // popis cinosti

		no = conp.total_particles();                             // pocet castic v kontejneru
		pst = try_add_2(id, nx, ny, nz, r, conp, theta, alfa, beta, B, fc);		  	// spocte pravdepodobnost se kterou dojde k operaci ADD
																				// pst = 0.5;
//		std::cout << pst << " ";
		pst = pst*(zet / (no + 1));
//		std::cout << pst << " ";

		// rn = rnd();							// generate random number between 0 and 1 
		if (rn2 < pst) {
			conp.put(id, nx, ny, nz, r);		// pridani castice
//			std::cout << "YES \n";
			if (id > npart) {				// pouzilo-li se ID o jedna vyssi nez maximalni ID v kontejneru a pridani
				npart = npart + 1;			// castice se uskutecni, pak je nutne navysit toto maximalni ID o 1
			}
			else {							// pouzilo-li se ID z fid, musi ve fid smazat
				fid.erase(fid.begin());
			}
			rn_add++;						// counter of added particles
		}
//		else { std::cout << "NO \n"; }
		//std::cout << npart << "\n";  /////////////////////////////////////////////////////////////////////////////////
	}     // BIRTH

	if (((const1) < rn1) && (rn1 <= (const2))) {
		del = uniform_int(1, conp.total_particles());	// vybere castici (resp. jeji poradi v containeru) ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
//		std::cout << "DELETE " << del << " " << conp.p[ijk][4 * q] << " " << conp.p[ijk][4 * q + 1] << " " << conp.p[ijk][4 * q + 2] << " " << conp.p[ijk][4 * q + 3] << "\n"; // popis cinosti
																																		   // nyni je potreba tuto castici odstranit z datovych struktur containeru, nestane se container nestabilni???
																																		   // pred vymazanim je treba si zapsat jeji sousedy:

		no = conp.total_particles();                     // pocet castic v kontejneru
		pst = try_delete_2(ijk, q, conp, theta, alfa, beta, B, fc);		// spocte pravdepodobnost se kterou dojde k operaci DELETE
																		// pst = 0.5;
//		std::cout << pst << " ";
		pst = pst*(no / zet);
//		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {
			erase(ijk, q, &fid, &conp);						// smazani castice - uvolni ID do fid k opetovnemu pouziti
//			std::cout << "YES \n";
			rn_del++;										// counter of deleted particles
		}
//		else { std::cout << "NO \n"; }
	}     // DEATH

	if (((const2) < rn1) && (rn1 <= (const3))) {
		del = uniform_int(1, conp.total_particles());	// vybere castici ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; // urci jeji souradnice
		r = conp.p[ijk][4 * q + 3];
		// vygeneruje novou castici z normalniho rozdeleni:
		nx = normal(x, sigma); ny = normal(y, sigma); nz = normal(z, sigma);   // coordinates of new particle
																			   //std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
																			   // these coordinates could be outside the window [0,1]^3 and in this case it is necessary to remap it into window:
		nx = nx - step_int(nx); ny = ny - step_int(ny); nz = nz - step_int(nz);
//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti

		pst = try_move_2(ijk, q, nx, ny, nz, conp, theta, alfa, beta, B, fc);		// urci pravdepodobnost se kterou dojde k operaci MOVE
																					// pst = 0.5;
//		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {
			erase(ijk, q, &conp);					// smazani castice
			conp.put(del, nx, ny, nz, r);				// pridani nove castice (ID zustava zachovano)

//			std::cout << "YES \n";
			rn_mov++;								// counter of moved particles
		}
//		else { std::cout << "NO \n"; }
	}     // MOVE 


	if (((const3) < rn1)) {
		del = uniform_int(1, conp.total_particles());	// vybere castici ke zmene (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
														//x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; // urci jeji souradnice
		r = conp.p[ijk][4 * q + 3];	// urci jeji polomer

		// vygeneruje novy polomer (3 moznosti):
		// 1) z useknuteho normalniho rozdeleni:
		//do {
		nr = normal(r, 0.0001);	// radius of new particle
		//} while (nr < 0);
		// 2) z exponencialniho rozdeleni
		//nr = exponential(50);			// double lambda = 50;
		// 3) z trohuhelnikoveho rozdeleni
		//nr = triangle(0.005, 0.03, 0.02);	// a = 0.005, b = 0.03, c = 0.02
		// pozn.: exponencialni rozdeleni dava nezaporne hodnoty, trojuhelnikove navic i omezene; zato normalni rozdeleni je na cele realne ose, pro 
		//	zjednoduseni ale uvazujme pouze jeho kladnou cast
		

//		std::cout << "CHANGE of RADIUS " << del << " " << r << " --> " << nr << "\n"; // popis cinosti
		if (nr < 0) {  } // not allowing negative radius
		else {
			pst = try_change_rad(ijk, q, nr, conp, theta, alfa, beta, B);		// urci pravdepodobnost se kterou dojde k operaci MOVE
																				// pst = 0.5;
//			std::cout << pst << " ";

																				// rn = rnd();  // generate random number between 0 and 1 
			if (rn2 <= pst) {
				conp.p[ijk][4 * q + 3] = nr;

//				std::cout << "YES \n";
				rn_chr++;								// counter of particles with changed radii
			}
//			else { std::cout << "NO \n"; }
		}

	}	  // CHANGE of RADIUS
}*/
