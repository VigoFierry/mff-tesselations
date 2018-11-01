
#include "Header.h"

// fcs_try version Laguerre (10.12.2017):
//		pridany varianty fci pro vypocet Laguerrova diagramu (navic fc try_change_rad); feasibilitu nutne resit pro zadanou skupinu bunek, zjednoduseni jako
//		u Voronoie neni mozne; upraven vypocet tezist (barycenter) a argumenty fce V2

// fcs_try version Fease (26.11.2017) - upusteno od: 
//              ve funkcich try_add, try_del jsou pozmeneny podminky pro overeni pripustnosti 
//              (motivace: pri pridani castice by nemela byt narusena horni mez, pri odebrani naopak dolni)
//              obdobna uprava pro fci try_move by mela nasledovat, zatim ponechana beze zmeny

// fcs_try version Addel (23.10.2017): 
//		ve funkcich try_add, try_del jsou od sebe oddeleny casti overujici feasibilitu a pocitajici strukturu, ci nektere casti energie 
//		(motivace: ve vetsine pripadu je pridani/odebrani bodu zastaveno pripustnosti)
//		fce try_move je ponechana, nebot je malokdy ukoncena nepripustnosti konfigurace

using namespace voro; // kvuli zavedeni trid (napr kvuli pouziti slova pre_container)
// using namespace std; - a nemusel bych psat vsude std::

					  
// fce typu TRY

// fce try_add "zkusi" pridat novou castici, tj. spocte a vrati hodnotu psti s kterou se bude pridavat nova castice v kroku ADD algoritu BDM
// potrebne promenne:	- id, nx, ny, nz (udaje o nove castici)
//						- con (container castic)
//						- alfa, beta, B, theta (parametry obecneho nastaveni)

// 18 - 380
double try_add(int id, double &nx, double &ny, double &nz, container &con, const double &theta, const double &alfa, const double &beta, const double &B)
{
	// [in]		id;x,y,z		the ID number and coordinates of the particle.
	// [in]		con				the container with stored particles.
	// [in]		theta			parameter
	// [in]		alfa,beta,B		hardcore parameters.

	unsigned int i, j, k, k1, k2, l, fng;
	int ijk, q, ijk_add, q_add, ijk_i, q_i, ijk_j, q_j;
	double part1, part2, part3, part4, part5;
	double xn, yn, zn, xnn, ynn, znn, vol, voll, h_max2, dist;
	bool in = false;

	std::vector<int> sr, tr;		// vectors storing ijk,q information about secondary, terciary particles
	std::vector<int> ntr;			// vector containing numbers of terciary particles for all secondary particles
	std::vector<int> sio, tio;		// information in/out for secondary, terciary particles
	std::vector<double> sap, tap;	// actual positions of secondary, terciary particles which are "out"
	std::vector<int> neigh, neighi;	// vector containing neighbor's ID
	std::vector<int> vert, verti;
	voronoicell_neighbor c, d;		// bunka s informaci o sousedech
	voronoicell di;

								//double alfa, beta; // asi je treba predat teto funkci parametricky
								//alfa = 0.02;
								//beta = 0.2;

	con.put(id, nx, ny, nz);   // pridani nove castice do containeru  
	find_pos(ijk_add, q_add, id, &con); // najde pozici nove castice 
										//
	// lepsi overit feasibility spolecne s vypoctem part1 (vypocet ma spolecne prvky)
	// Addel: pokus zda oddeleni techto casti nezvysi rychlost simulace (ve vetsine pripadu je pridani bodu zastaveno pripustnosti)
	
	con.compute_cell(c, ijk_add, q_add);	// {1,F} compute cell for new added particle

	c.neighbors(neigh);						// {1,F} compute its neighbors			
	c.face_vertices(vert);					// {1PER} computes list of vertices for each face of the cell
	vol = c.volume();
	fng = 0;								// {1PER}

	
	// A  - overeni pripustnosti
	/*
	// jiny pristup feasibility: - zde neni moc vyhodny, protoze bych potreboval vsechny zmenene castice, tj. pr, sr
	cells = sr;
	cells.push_back(ijk_add);
	cells.push_back(q_add);

	std::cout << "NEW Feasibility: ";
	if (feasibility(con, cells, alfa, beta, B)) { std::cout << "OK. \n"; }
	else { con.put(del, x_del, y_del, z_del); return 0; }
	//-----------------------------------------------------------------------------------*/

	// number of faces
	if (c.number_of_faces() > 12) { erase(ijk_add, q_add, &con); return 0; };


	for (i = 0; i < neigh.size(); i++) {	// {1,F} loop over the neighbors
		find_pos(ijk, q, neigh[i], &con);   // {1,F} find position of the i-th neighbor
											// std::cout << con.p[ijk][3 * q] << " " << con.p[ijk][3 * q + 1] << " " << con.p[ijk][3 * q + 2] << " ; "; /////////////////

		face_dist(fng, vert, nx, ny, nz, xn, yn, zn, c);		// {1PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
		xn = nx + 2 * xn; yn = ny + 2 * yn; zn = nz + 2 * zn;	// {1PER} urcim skutecne souradnice tohoto souseda za stenou

		// std::cout << xn << " " << yn << " " << zn << " || "; //////////////////////////////////////////////////////////
		// - feasibility musi kalkulovat se "skutecnymi" souradnicemi !!!!!!!!!!!!!
		dist = point_dist(nx, ny, nz, xn, yn, zn) / 2;

		//if (dist > beta) { std::cout << "A: Tesselation is NOT feasible (" << point_dist(nx, ny, nz, xn, yn, zn) << " - upper bound).\n"; erase(ijk_add, q_add, &con); return 0; }; 
		//if (dist < alfa) { std::cout << "A: Tesselation is NOT feasible (" << point_dist(nx, ny, nz, xn, yn, zn) << " - lower bound).\n"; erase(ijk_add, q_add, &con); return 0; };

		if (dist > beta) { erase(ijk_add, q_add, &con); return 0; };
		if (dist < alfa) { erase(ijk_add, q_add, &con); return 0; };

		if (B > 0) {  // B=0 znamena bez omezeni tvaru 
			con.compute_cell(d, ijk, q);        // {1} compute the cell of this neighbor

			// nova bunka:
			if (d.number_of_faces() > 12) { erase(ijk_add, q_add, &con); return 0; };

			//if ((pow(dist, 3) / vol) > B) { std::cout << "A: Tesselation is NOT feasible (shape constraint - " << (pow(dist, 3) / vol) << " vs. " << B << ").\n"; erase(ijk_add, q_add, &con); return 0; };
			if ((pow(dist, 3) / vol) > B) { erase(ijk_add, q_add, &con); return 0; };

			// pridanim bodu se zmensi objemy jeho sousedu, je proto nutne overit podminku pro vsechny bunky, jejichz objem se zmenil:
			d.neighbors(neighi);						// a jeho sousedy
			d.face_vertices(verti);						// computes list of vertices for each face of the cell
			voll = d.volume();
			j = neighi.size();
			h_max2 = pow(h_maximum(j, verti, d, xn, yn, zn), 3);

			//if (h_max2 > B*voll) { std::cout << "A: Tesselation is NOT feasible (shape constraint - " << h_max2 / voll << " vs. " << B << ").\n"; erase(ijk_add, q_add, &con); return 0; }; // double sqrt(double x);			
			if (h_max2 > B*voll) { erase(ijk_add, q_add, &con); return 0; }; 		

			//if (h_max2 > 4 * B*voll) { /*std::cout << "A: Tesselation is NOT feasible (shape constraint - " << h_max2 << " vs. " << 4 * B*voll << ").\n";*/ erase(ijk_add, q_add, &con); return 0; }; // double sqrt(double x);			
		}
		fng = fng + vert[fng] + 1;			// {1PER} set actual position in vector of face vertices
	}

	k1 = 0;
	// spocti pravdepodobnost s kterou pridame castici; nejdrive spoctem vnitrek exponenciely po castech:
	//		part1, part2, part3  se tykaji bunek po pridani nove castice
	part1 = 0; // V2 pro nove pridanou bunku a jeji primarni sousedy (secondary particles)
	fng = 0;

	// B  - ulozeni struktury sekundarnich castic a vypocet part 1
	for (i = 0; i < neigh.size(); i++) {	// {1,F} loop over the neighbors
		find_pos(ijk, q, neigh[i], &con);   // {1,F} find position of the i-th neighbor
											// std::cout << con.p[ijk][3 * q] << " " << con.p[ijk][3 * q + 1] << " " << con.p[ijk][3 * q + 2] << " ; "; /////////////////

		face_dist(fng, vert, nx, ny, nz, xn, yn, zn, c);		// {1PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
		xn = nx + 2 * xn; yn = ny + 2 * yn; zn = nz + 2 * zn;	// {1PER} urcim skutecne souradnice tohoto souseda za stenou

		// std::cout << xn << " " << yn << " " << zn << " || "; //////////////////////////////////////////////////////////
		// - kdyby byla podminka if uz zde, a kdyby se pri vysledku IN zpresnila hodnota xn == con.p[ijk][3 * q], ...
		//        byl by test feasibility presnejsi  
		
		sr.push_back(ijk); sr.push_back(q); // add this neighbor to the second particle vector
		con.compute_cell(d, ijk, q);        // {1} compute the cell of this neighbor

		// TEST		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) { std::cout << "IN 1 " << xn << " " << con.p[ijk][3 * q] << " " << yn << " " << con.p[ijk][3 * q + 1] << " " << zn << " " << con.p[ijk][3 * q + 2] << "\n"; }
		// TEST		else { std::cout << "OUT 1 " << xn << " " << con.p[ijk][3 * q] << " " << yn << " " << con.p[ijk][3 * q + 1] << " " << zn << " " << con.p[ijk][3 * q + 2] << "\n"; }
		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		 // {1PER} je-li pozice souseda v okne (tj. nedoslo k preklopeni)

			// std::cout << "IN \n"; ///////////////////////////////////////////////////////////////////////////////////
			part1 = part1 + V2(c, d, nx, ny, nz, xn, yn, zn);    // {1} add the value of V2 function to the first sum
			sio.push_back(0);					// neighbor is "in"
		}
		else {
			part1 = part1 + V2(c, d, nx, ny, nz, xn, yn, zn);
			k1++;
			sio.push_back(k1);					// neighbor is k-th "out"
			sap.push_back(xn); sap.push_back(yn); sap.push_back(zn);  // storing real position of this neighbor
		}
		// std::cout << V2(theta, c, d) << " \n"; ///////////////////////////////////////////////////////////////////////
		fng = fng + vert[fng] + 1;			// {1PER} set actual position in vector of face vertices
	}
	/*std::cout << " \n"; /////////////////////////////////////////////////////////////////////////////////////////////////
	for (i = 0; i < sr.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	std::cout << sr[2 * i] << " " << sr[2 * i + 1] << " ; ";
	} std::cout << " \n";
	for (i = 0; i < sio.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	std::cout << sio[i] << " ";
	} std::cout << " \n";
	for (i = 0; i < sap.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	std::cout << sap[i] << " " << sap[i + 1] << " " << sap[i + 2] << " ; ";
	} std::cout << " \n"; */

//	std::cout << "part1: " << part1 << " \n"; ///////////////////////////////////////////////////////////////////////////
											  // std::cout << sr.size() << " sr \n";  /////////////////////////////////////////////////////////////////////////////
											  // for (i = 0; i < sr.size() / 2; i++) { std::cout << sr[2 * i] << " " << sr[2 * i + 1] << " \n"; } //////////////////

	k2 = 0; l = 0;
	part2 = 0; // V2 pro primarni sousedy (po pridani)
	part3 = 0; // V2 pro primarni a sekundarni sousedy (po pridani)
	for (i = 0; i < sr.size() / 2; i++) {                           // {2,3} loop over secondary particles
		ijk_i = sr[2 * i]; q_i = sr[2 * i + 1];

		if (sio[i] == 0) {
			xnn = con.p[ijk_i][3 * q_i];							// {2,3PER} computes coordinates of i-th particle
			ynn = con.p[ijk_i][3 * q_i + 1];
			znn = con.p[ijk_i][3 * q_i + 2];
		}
		else {
			xnn = sap[3 * (sio[i] - 1)];
			ynn = sap[3 * (sio[i] - 1) + 1];
			znn = sap[3 * (sio[i] - 1) + 2];
		}

		con.compute_cell(c, ijk_i, q_i);				// {2,3} compute the cell of the i-th secundary particle

		for (j = i + 1; j < sr.size() / 2; j++) {					// {2} loop over secondary particles with "higher order" (to avoid double counting)  	
			ijk_j = sr[2 * j]; q_j = sr[2 * j + 1];
			if (are_neighbors(c, ijk_j, q_j, &con)) {	  // if i and j are neighbors... (i.e. if c_i ~ c_j)
				con.compute_cell(d, ijk_j, q_j);		  // compute cell j

				if (sio[j] == 0) {
					xn = con.p[ijk_j][3 * q_j];							// {2,3PER} computes coordinates of i-th particle
					yn = con.p[ijk_j][3 * q_j + 1];
					zn = con.p[ijk_j][3 * q_j + 2];
				}
				else {
					xn = sap[3 * (sio[j] - 1)];
					yn = sap[3 * (sio[j] - 1) + 1];
					zn = sap[3 * (sio[j] - 1) + 2];
				}

				part2 = part2 + V2(c, d, xnn, ynn, znn, xn, yn, zn);
			}
		}
		
		c.neighbors(neigh);											// {3} compute its neighbors
		c.face_vertices(vert);										// {3PER} computes list of vertices for each face of the cell
		fng = 0;													// {3PER}

		for (j = 0; j < neigh.size(); j++) {						// {3} loop over these neighbors

			if (terciary(neigh[j], id, sr, &con)) {			          // {3} je-li soused "terciarni castice", tj je sousedem souseda pridane castice a zaroven neni sousedem pridane castice
				find_pos(ijk, q, neigh[j], &con);					    // najdi pozici souseda
				tr.push_back(ijk); tr.push_back(q);					    // uloz souseda mezi terciarni 

				// !!!!! prvni tri souradnice dane fci face_dist musi odpovidat stredu bunky, tj. musi byt z [0,1]^3
				face_dist(fng, vert, con.p[ijk_i][3 * q_i], con.p[ijk_i][3 * q_i + 1], con.p[ijk_i][3 * q_i + 2], xn, yn, zn, c);		// {3PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
				xn = xnn + 2 * xn; yn = ynn + 2 * yn; zn = znn + 2 * zn;// {3PER} urcim skutecne souradnice tohoto souseda za stenou

																		// na NIC			if (xn == con.p[ijk][3 * q] && yn == con.p[ijk][3 * q + 1] && zn == con.p[ijk][3 * q + 2]) { // {3PER} nedoslo k preklopeni (terciarni castice je "in")
				if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {  // {3PER} nedoslo k preklopeni (terciarni castice je "in")
					tio.push_back(0);
					con.compute_cell(d, ijk, q);
					part3 = part3 + V2(c, d, xnn, ynn, znn, xn, yn, zn);
				}
				else {								// doslo k preklopeni (terciarni catice je "out")
					k2++; tio.push_back(k2);
					tap.push_back(xn); tap.push_back(yn); tap.push_back(zn); // storing real position of this neighbor

					con.compute_cell(d, ijk, q);
					part3 = part3 + V2(c, d, xnn, ynn, znn, xn, yn, zn);
						/*part3 = part3 + V2(c, d, false, con.p[ijk][3 * q], con.p[ijk][3 * q + 1], con.p[ijk][3 * q + 2],
							xn, yn, zn);*/
						// con.p[sr[2 * i]][3 * sr[2 * i + 1]], con.p[sr[2 * i]][3 * sr[2 * i + 1] + 1], con.p[sr[2 * i]][3 * sr[2 * i + 1] + 2], xn, yn, zn);
						
				} // END if..else (preklopeni)
				l++;
			} // END if (terciary)
			fng = fng + vert[fng] + 1;			// {2,3PER} set actual position in vector of face vertices (posun na dalsi stenu)
		} // END for (j; cyklus over neighbors of secondary particle i)
		ntr.push_back(l);									// uloz pocet terciarnich castic souseda i
	} // END for (i; cyklus over sr)
//	std::cout << "part2: " << part2 << " \n";  ///////////////////////////////////////////////////////////////////// 

											   // std::cout << terciary(5, id, sr, &con) << " \n"; ///////////////////////////////////////////////////////////////
//	std::cout << "part3: " << part3 << " \n";  /////////////////////////////////////////////////////////////////////

	/*std::cout << tr.size() << " tr \n";  ///////////////////////////////////////////////////////////////////////////
	 std::cout << ntr.size() << " ntr \n"; //////////////////////////////////////////////////////////////////////////
	 for (i = 0; i < ntr.size(); i++) {
	 std::cout << ntr[i] << " ";
	 } std::cout << " \n";  ///////////////////////////////////////////////////////////////////////////////////////////////
	 for (i = 0; i < tr.size() / 2; i++) {
	 std::cout << tr[2 * i] << " " << tr[2 * i + 1] << " ; ";
	} std::cout << " \n";  /////////////////////////////////////////////////////////////////////////////////////////
	 for (i = 0; i < tio.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
											   std::cout << tio[i] << " ";
											   } std::cout << " \n";
											   for (i = 0; i < tap.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
											   std::cout << tap[i] << " " << tap[i + 1] << " " << tap[i + 2] << " ; ";
											   } std::cout << " \n"; */
											   // std::cout << con.total_particles() << " \n"; ///////////////////////////////////////////////////////////////////

											   // nyni jsme spocetli cast vnitrku exponenciely tykajici se stavu po pridani bunky;
											   // dale zbyva cast tykajici se stavu pred pridanim 
											   // - odstranim novou castici a vyuziju jiz spoctenych sekundarnich a terciarnich castic
											   // - nevyhodou je ze musim znovu cyklit, toho by se slo zbavit tim, ze bych mel ulozene dva containery - pro stavy
											   //	 pred a po  (je dulezitejsi casova uspora oproti pametove?)

	erase(ijk_add, q_add, &con); // smaz novou castici at muzeme spocitat stav pred pridanim

								 // std::cout << con.total_particles() << " \n";  //////////////////////////////////////////////////////////////////

								 // part4, part5  se tykaji bunek pred pridanim nove castice
	part4 = 0; // V2 pro primarni sousedy (pred pridanim)  // cast je totozna s part2
	part5 = 0; // V2 pro primarni a sekundarni sousedy (pred pridanim)
	for (i = 0; i < sr.size() / 2; i++) {							// {4,5} loop over secondary particles
		ijk_i = sr[2 * i]; q_i = sr[2 * i + 1];
		con.compute_cell(c, ijk_i, q_i);							// {4,5} compute its cells
																	// c.neighbors(neigh);										// {4,5} compute its neighbors
		if (sio[i] == 0) {
			xn = con.p[ijk_i][3 * q_i];								// {F} computes coordinates of i-th particle
			yn = con.p[ijk_i][3 * q_i + 1];
			zn = con.p[ijk_i][3 * q_i + 2];
		}
		else {
			xn = sap[3 * (sio[i] - 1)];
			yn = sap[3 * (sio[i] - 1) + 1];
			zn = sap[3 * (sio[i] - 1) + 2];
		}

		for (j = i + 1; j < sr.size() / 2; j++) {					// {4} loop over secondary particles with "higher order" (to avoid double counting)  	
			ijk_j = sr[2 * j]; q_j = sr[2 * j + 1];
			if (are_neighbors(c, ijk_j, q_j, &con)) {	  // if i and j are neighbors... (i.e. if c_i ~ c_j)
				con.compute_cell(d, ijk_j, q_j);

				if (sio[j] == 0) {
					xnn = con.p[ijk_j][3 * q_j];								// {F} computes coordinates of i-th particle
					ynn = con.p[ijk_j][3 * q_j + 1];
					znn = con.p[ijk_j][3 * q_j + 2];
				}
				else {
					xnn = sap[3 * (sio[j] - 1)];
					ynn = sap[3 * (sio[j] - 1) + 1];
					znn = sap[3 * (sio[j] - 1) + 2];
				}

				part4 = part4 + V2(c, d, xn, yn, zn, xnn, ynn, znn);
			}
		}

		if (i == 0) {
			for (j = 0; j < ntr[0]; j++) {							// {5} i=0; loop over its terciary particles
				ijk_j = tr[2 * j]; q_j = tr[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		  // compute cell of this terciary particle

				if (tio[j] == 0) {
					xnn = con.p[ijk_j][3 * q_j];								// {F} computes coordinates of i-th particle
					ynn = con.p[ijk_j][3 * q_j + 1];
					znn = con.p[ijk_j][3 * q_j + 2];
				}
				else {
					xnn = tap[3 * (tio[j] - 1)];
					ynn = tap[3 * (tio[j] - 1) + 1];
					znn = tap[3 * (tio[j] - 1) + 2];
				}

				part5 = part5 + V2(c, d, xn, yn, zn, xnn, ynn, znn);
			}
		}
		else {
			for (j = 0; j < ntr[i] - ntr[i - 1]; j++) {				// {5} i>0; loop over its terciary particles
				k = ntr[i - 1] + j;
				ijk_j = tr[2 * k]; q_j = tr[2 * k + 1];
				con.compute_cell(d, ijk_j, q_j);

				if (tio[k] == 0) {
					xnn = con.p[ijk_j][3 * q_j];								// {F} computes coordinates of i-th particle
					ynn = con.p[ijk_j][3 * q_j + 1];
					znn = con.p[ijk_j][3 * q_j + 2];
				}
				else {
					xnn = tap[3 * (tio[k] - 1)];
					ynn = tap[3 * (tio[k] - 1) + 1];
					znn = tap[3 * (tio[k] - 1) + 2];
				}

				part5 = part5 + V2(c, d, xn, yn, zn, xnn, ynn, znn);

			
				
			}
		} // END if..else (i==0 vs i>0)

	} // END for (i; cyklus over secondary particles)

//	std::cout << "part4: " << part4 << " \n";  ///////////////////////////////////////////////////////////////////// 
//	std::cout << "part5: " << part5 << " \n";  ///////////////////////////////////////////////////////////////////// 

//	std::cout << "A: part1: " << -part1 << " | part2: " << -part2 << " | part3: " << -part3 << " | part4: " << part4 << " | part5: " << part5 << " \n";

	return exp(theta*(-part1 - part2 - part3 + part4 + part5));
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
// _____________________________________________________________________________________________________________________________________________


// ---------------------------------------------------------------------------------------------------------------------------------------------
// _____________________________________________________________________________________________________________________________________________
// _____________________________________________________________________________________________________________________________________________



// 395 - 755

double try_delete(int &ijk_del, int &q_del, container &con, const double &theta, const double &alfa, const double &beta, const double &B)
{
	// [in]		ijk_del,q_del	position of the particle.
	// [in]		con				the container with stored particles.
	// [in]		theta			parameter
	// [in]		alfa,beta,B		hardcore parameters.

	unsigned int i, j, k, k1, k2, l, fng;
	int ijk, q, ijk_i, q_i, ijk_j, q_j, del;
	double part1, part2, part3, part4, part5;
	double x_del, y_del, z_del, x, y, z, xn, yn, zn;

	std::vector<int> pr, sr, tr; // misto particle_order class snadnejsi pouzivat vektory pro ukladani ijk,q castic
	std::vector<int> sio, tio;
	std::vector<double> sap, tap;
	std::vector<int> ntr, cells;
	std::vector<int> neigh;  // vektor pro prirazeni ID sousedu
	std::vector<int> vert;
	voronoicell_neighbor c_del, c, d;  // bunka s informaci o sousedech

									   //double alfa, beta; // asi je treba predat teto funkci parametricky
									   //alfa = 0.05;
									   //beta = 0.5;

									   // priprava: najdu souradnice a ID mazane castice, spoctu jeji bunku, pak jeji sousedy a vrcholy sten
									   //		loop pres sousedy ulozi do sr pozice ijk,q techto sekundarnich castic, face_dist spocte skutecne polohy 
									   //		techto sousedu, do sio se ulozi informace o preklopeni a do sap pripadne skutecna poloha sekundarni castice

	x_del = con.p[ijk_del][3 * q_del];      // coordinates of deleted particle
	y_del = con.p[ijk_del][3 * q_del + 1];
	z_del = con.p[ijk_del][3 * q_del + 2];
	del = con.id[ijk_del][q_del];			// ID of deleted particle

	con.compute_cell(c_del, ijk_del, q_del);// {} compute the cell of deleted particle
	c_del.neighbors(neigh);					// {} compute its neighbors (in following "primary")
	c_del.face_vertices(vert);				// {PER} computes list of vertices for each face of the cell
	fng = 0;								// {PER}

											//std::cout << del << " " << ijk_del << "\n";  /////////////////////////////////////////////////////////////////////////
											//std::cout << "uvod /n";  /////////////////////////////////////////////////////////////////////////////////////////////
											// nejdrive se spocetli udaje o mazane castici, jeji bunka, sousede, vrcholy sten; pak se castice smaze, a az pak se 
											//    ulozi vektory sr,sio,sap (udaje o castici c_del zustaly ulozeny, takze nevadi ze tato bunka uz nebude existovat)
											//    (nutne je toto poradi, protoze fce erase meni hodnoty q, tj pozice v boxu)

											// smazani castice a vypocet stavu po smazani (mezi prvnimi nas zajima feasibilita, bez ni nema smysl pokracovat)
	erase(ijk_del, q_del, &con); // smaz castici at muzeme spocitat stav po smazani, smazani zmeni strukturu

								 //std::cout << "smazano \n"; ///////////////////////////////////////////////////////////////////////////////////////////

	k1 = 0;
	for (i = 0; i < neigh.size(); i++) {	// {} loop over the neighbors of deleted particle
		find_pos(ijk, q, neigh[i], &con);   // {} find position of the i-th primary neighbor after deleting 
		sr.push_back(ijk); sr.push_back(q); // add this primary neighbor to the vector of secondary particles

		face_dist(fng, vert, x_del, y_del, z_del, xn, yn, zn, c_del);	// {PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
		xn = x_del + 2 * xn; yn = y_del + 2 * yn; zn = z_del + 2 * zn;	// {PER} urcim skutecne souradnice tohoto souseda za stenou

																		//		if (xn == con.p[ijk][3 * q] && yn == con.p[ijk][3 * q + 1] && zn == con.p[ijk][3 * q + 2]) { // {PER} nedoslo k preklopeni
		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {PER} nedoslo k preklopeni
			sio.push_back(0);					// neighbor is "in"
		}
		else {
			k1++;
			sio.push_back(k1);					// neighbor is k-th "out"
			sap.push_back(xn); sap.push_back(yn); sap.push_back(zn);  // storing real position of this neighbor
		}

		fng = fng + vert[fng] + 1;
	} // END for (i; cyklus over neighbors of deleted particle)

	  /* for (i = 0; i < neigh.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << neigh[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sr.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << sr[2 * i] << " " << sr[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < sio.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sio[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sap.size()/3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sap[i] << " " << sap[i + 1] << " " << sap[i + 2] << " ; ";
	  } std::cout << " \n"; */
	  //std::cout << "priprava /n"; //////////////////////////////////////////////////////////////////////////////////////////

	// Addel: striktni oddeleni overeni pripustnosti od ostatnich ukonu, pripustnost je prvorada

	
	// jiny pristup feasibility:
	cells = sr;

//	std::cout << "NEW Feasibility: ";
	if (feasibility(con, cells, alfa, beta, B)) { /*std::cout << "OK. \n";*/ }
	else { con.put(del, x_del, y_del, z_del); return 0; }
	//-----------------------------------------------------------------------------------

	

	// B  - vypocet struktury terciarnich castic a part 1, part 2 
	part1 = 0;
	part2 = 0;
	l = 0;
	k2 = 0;
	for (i = 0; i < sr.size() / 2; i++) {                           // {1,2,F} loop over secondary particles
		ijk_i = sr[2 * i]; q_i = sr[2 * i + 1];
		con.compute_cell(c, ijk_i, q_i);							// {1,2,F} compute cell of i-th secondary particle

		if (sio[i] == 0) {
			x = con.p[ijk_i][3 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][3 * q_i + 1];
			z = con.p[ijk_i][3 * q_i + 2];
		}
		else {
			x = sap[3 * (sio[i] - 1)];
			y = sap[3 * (sio[i] - 1) + 1];
			z = sap[3 * (sio[i] - 1) + 2];
		}

		for (j = i + 1; j < sr.size() / 2; j++) {                   // {1,F} loop over secondary particles with "higher order" (to avoid double counting) 
			if (are_neighbors(c, sr[2 * j], sr[2 * j + 1], &con)) {	// {1,F} sousedi-li bunky i a j ... (i.e. if c_i ~ c_j)
				ijk_j = sr[2 * j]; q_j = sr[2 * j + 1];
				// std::cout << "are \n"; ///////////////////////////////////////////////////////////////////////////////
				if (sio[j] == 0) {
					xn = con.p[ijk_j][3 * q_j];				// {F} computes coordinates of j-th neighbor
					yn = con.p[ijk_j][3 * q_j + 1];
					zn = con.p[ijk_j][3 * q_j + 2];
				}
				else {
					xn = sap[3 * (sio[j] - 1)];
					yn = sap[3 * (sio[j] - 1) + 1];
					zn = sap[3 * (sio[j] - 1) + 2];
				}
			
				con.compute_cell(d, ijk_j, q_j);					// {1} compute cell of the j-th secondary particle

				part1 = part1 + V2(c, d, x, y, z, xn, yn, zn);
				
			} // END if (are_neighbors)
		} // END for (j; secondary particles with higher order - j>i)
		c.neighbors(neigh);						// {2} compute its neighbors ("secondary")
		c.face_vertices(vert);					// {2PER} computes list of vertices for each face of the cell
		fng = 0;								// {2PER}


		for (j = 0; j < neigh.size(); j++) {				           // {2} loop over terciary particles (secondary neighbors)
			if (terciary(neigh[j], sr, &con)) {							 // je-li castice terciarni, tj soused sekundarni castice, ktery zaroven neni sekundarni castici
				find_pos(ijk, q, neigh[j], &con);			             // find its position
				tr.push_back(ijk); tr.push_back(q);			             // store this secondary neighbor to the vector of terciary particles 

																		 // pro volani fce face_dist musim pouzit stred bunky, tj souradnice ktere jsou v okne !!! a nedbat na to, zda jsou spravne
				face_dist(fng, vert, con.p[ijk_i][3 * q_i], con.p[ijk_i][3 * q_i + 1], con.p[ijk_i][3 * q_i + 2], xn, yn, zn, c);		// {3PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
																																		// dale jiz ale mohu pouzit skutecne souradnice bodu bunky (x,y,z):
				xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;// {3PER} urcim skutecne souradnice tohoto souseda za stenou

																  //					if (xn == con.p[ijk][3 * q] && yn == con.p[ijk][3 * q + 1] && zn == con.p[ijk][3 * q + 2]) { // {2PER} nedoslo k preklopeni terciarni castice (terciarni castice je "in")
				if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {2PER} nedoslo k preklopeni terciarni castice (terciarni castice je "in")
					tio.push_back(0);
					con.compute_cell(d, ijk, q);
					part2 = part2 + V2(c, d, x, y, z, xn, yn, zn);

					
				}
				else {								// terciarni castice je preklopena
					k2++; tio.push_back(k2);
					tap.push_back(xn); tap.push_back(yn); tap.push_back(zn); // storing real position of this neighbor
					
					con.compute_cell(d, ijk, q);
					part2 = part2 + V2(c, d, x, y, z, xn, yn, zn);
					
				} // END if..else (terciary particle is IN vs OUT)
				l++;
			} // END if (terciary)
			fng = fng + vert[fng] + 1;
		}  // END for (j; loop over neighbors of i-th secondary particle)
		ntr.push_back(l);									         // {2} uloz pocet sekundarnich sousedu souseda i
	} // END for (i; loop over secondary particles)

//	std::cout << "part1: " << part1 << " \n";  //////////////////////////////////////////////////////////////////////
//	std::cout << "part2: " << part2 << " \n";  //////////////////////////////////////////////////////////////////////

											   /* for (i = 0; i < ntr.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
											   std::cout << ntr[i] << " ";
											   } std::cout << " \n";
											   for (i = 0; i < tr.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
											   std::cout << tr[2 * i] << " " << tr[2 * i + 1] << " ; ";
											   } std::cout << " \n";
											   for (i = 0; i < tio.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
											   std::cout << tio[i] << " ";
											   } std::cout << " \n";
											   for (i = 0; i < tap.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
											   std::cout << tap[i] << " " << tap[i + 1] << " " << tap[i + 2] << " ; ";
											   } std::cout << " \n"; */


	con.put(del, x_del, y_del, z_del);
	// castice se prida sice asi do stejneho boxu, ale rozhodne ne na sve puvodni misto (nybrz na konec) !!!!!!!!
	//    !!!!!!! je proto potreba urcit opravdovou polohu teto castice !!!!!!!!!!!!!
	find_pos(ijk_del, q_del, del, &con);

	part3 = 0;
	for (i = 0; i < sr.size() / 2; i++) {				// {3} loop over the secondary particles
		ijk_i = sr[2 * i]; q_i = sr[2 * i + 1];
		con.compute_cell(d, ijk_i, q_i);  // {3} compute the cell of this primary neighbor

		if (sio[i] == 0) {
			x = con.p[ijk_i][3 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][3 * q_i + 1];
			z = con.p[ijk_i][3 * q_i + 2];
		}
		else {
			x = sap[3 * (sio[i] - 1)];
			y = sap[3 * (sio[i] - 1) + 1];
			z = sap[3 * (sio[i] - 1) + 2];
		}

		part3 = part3 + V2(c_del, d, x_del, y_del, z_del, x, y, z);


	}
//	std::cout << "part3: " << part3 << " \n";  /////////////////////////////////////////////////////////////////////

	part4 = 0;
	part5 = 0;
	for (i = 0; i < sr.size() / 2; i++) {                              // {4,5} loop over secondary particles    
		ijk_i = sr[2 * i]; q_i = sr[2 * i + 1];
		if (sio[i] == 0) {
			x = con.p[ijk_i][3 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][3 * q_i + 1];
			z = con.p[ijk_i][3 * q_i + 2];
		}
		else {
			x = sap[3 * (sio[i] - 1)];
			y = sap[3 * (sio[i] - 1) + 1];
			z = sap[3 * (sio[i] - 1) + 2];
		}
		con.compute_cell(c, ijk_i, q_i);				// {4,5} compute the cell of i-th primary neighbor   

		for (j = i + 1; j < sr.size() / 2; j++) {					// {4} loop over secondary particles with "higher order" (to avoid double counting)  	
			if (are_neighbors(c, sr[2 * j], sr[2 * j + 1], &con)) {	  // if i and j are neighbors... (i.e. if c_i ~ c_j)
				ijk_j = sr[2 * j]; q_j = sr[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);					  // compute cell j

				if (sio[j] == 0) {
					xn = con.p[ijk_j][3 * q_j];					// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][3 * q_j + 1];
					zn = con.p[ijk_j][3 * q_j + 2];
				}
				else {
					xn = sap[3 * (sio[j] - 1)];
					yn = sap[3 * (sio[j] - 1) + 1];
					zn = sap[3 * (sio[j] - 1) + 2];
				}

				part4 = part4 + V2(c, d, x, y, z, xn, yn, zn);
			}
		}

		if (i == 0) {
			for (j = 0; j < ntr[0]; j++) {							// {5} i=0; loop over its terciary particles
				ijk_j = tr[2 * j]; q_j = tr[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		  // compute cell of this terciary particle
				if (tio[j] == 0) {
					xn = con.p[ijk_j][3 * q_j];								// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][3 * q_j + 1];
					zn = con.p[ijk_j][3 * q_j + 2];
				}
				else {
					xn = tap[3 * (tio[j] - 1)];
					yn = tap[3 * (tio[j] - 1) + 1];
					zn = tap[3 * (tio[j] - 1) + 2];
				}

				part5 = part5 + V2(c, d, x, y, z, xn, yn, zn);
			}
		}
		else {
			for (j = 0; j < ntr[i] - ntr[i - 1]; j++) {				// {5} i>0; loop over its terciary particles
				k = ntr[i - 1] + j;
				ijk_j = tr[2 * k]; q_j = tr[2 * k + 1];
				con.compute_cell(d, ijk_j, q_j);

				if (tio[k] == 0) {
					xn = con.p[ijk_j][3 * q_j];								// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][3 * q_j + 1];
					zn = con.p[ijk_j][3 * q_j + 2];
				}
				else {
					xn = tap[3 * (tio[k] - 1)];
					yn = tap[3 * (tio[k] - 1) + 1];
					zn = tap[3 * (tio[k] - 1) + 2];
				}

				part5 = part5 + V2(c, d, x, y, z, xn, yn, zn);
			}
		} // END if..else (i==0 vs i>0)

	}  // END for (i; loop over secondary particles)

//	std::cout << "part4: " << part4 << " \n";  /////////////////////////////////////////////////////////////////////
//	std::cout << "part5: " << part5 << " \n";  /////////////////////////////////////////////////////////////////////

//	std::cout << "D: part1: " << -part1 << " | part2: " << -part2 << " | part3: " << part3 << " | part4: " << part4 << " | part5: " << part5 << " \n";

	return exp(theta*(-part1 - part2 + part3 + part4 + part5));
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
// _____________________________________________________________________________________________________________________________________________


// ---------------------------------------------------------------------------------------------------------------------------------------------
// _____________________________________________________________________________________________________________________________________________
// _____________________________________________________________________________________________________________________________________________



// 770 - 1227

double try_move(int &ijk_del, int &q_del, double &nx, double &ny, double &nz, container &con, const double &theta, const double &alfa, const double &beta, const double &B)
{
	// [in]		ijk_del,q_del	position of the particle to be deleted.
	// [in]		x,y,z			the ID number and coordinates of the particle to be added.  -- misto tohoto by stacilo predat hodnotu sigma
	// [in]		con				the container with stored particles.
	// [in]		theta			parameter
	// [in]		alfa,beta,B		hardcore parameters.

	unsigned int i, j, k, l, k1_add, k1_del, k2, fng;
	int ijk, q, ijk_add, q_add, del, ijk_i, q_i, ijk_j, q_j;
	double part1, part2, part3, part4, part5, part6;
	double x_del, y_del, z_del, x, y, z, xn, yn, zn;
	double vol;

	std::vector<int> sr_add, sr_del, sr, tr; // misto particle_order class snadnejsi pouzivat vektory pro ukladani ijk,q castic
	std::vector<int> sio_add, sio_del, tio;
	std::vector<double> sap_add, sap_del, tap;
	std::vector<int> ntr, cells;
	std::vector<int> neigh_del, neigh_add, neigh, neighi;  // vektor pro prirazeni ID sousedu
	std::vector<int> vert, verti;
	voronoicell_neighbor c_del, c_add, c, d;  // bunka s informaci o sousedech


	x_del = con.p[ijk_del][3 * q_del];      // coordinates of deleted particle
	y_del = con.p[ijk_del][3 * q_del + 1];
	z_del = con.p[ijk_del][3 * q_del + 2];

	del = con.id[ijk_del][q_del];			// ID of deleted particle

	con.compute_cell(c_del, ijk_del, q_del);	// compute the cell of deleted particle
	c_del.neighbors(neigh_del);					// compute its neighbors (in following "primary")
												//n_del = neigh_del.size();					// determine number of neighbors of deleted particle
	c_del.face_vertices(vert);					// {PER} computes list of vertices for each face of the cell
	fng = 0;

	erase(ijk_del, q_del, &con);				// delete old particle
												// std::cout << "del \n";     //////////////////////////////////////////////////////////////////////////////////////////////

	k1_del = 0;
	for (i = 0; i < neigh_del.size(); i++) {
		find_pos(ijk, q, neigh_del[i], &con);
		sr_del.push_back(ijk); sr_del.push_back(q);
		// std::cout << con.p[ijk][3 * q] << " " << con.p[ijk][3 * q + 1] << " " << con.p[ijk][3 * q + 2] << " ; "; /////////////////

		face_dist(fng, vert, x_del, y_del, z_del, xn, yn, zn, c_del);	// {PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
		xn = x_del + 2 * xn; yn = y_del + 2 * yn; zn = z_del + 2 * zn;	// {PER} urcim skutecne souradnice tohoto souseda za stenou

																		//		if (xn == con.p[ijk][3 * q] && yn == con.p[ijk][3 * q + 1] && zn == con.p[ijk][3 * q + 2]) { // {PER} nedoslo k preklopeni
		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {PER} nedoslo k preklopeni
			sio_del.push_back(0);					// neighbor is "in"
		}
		else {
			k1_del++;
			sio_del.push_back(k1_del);					// neighbor is k-th "out"
			sap_del.push_back(xn); sap_del.push_back(yn); sap_del.push_back(zn);  // storing real position of this neighbor
		}

		fng = fng + vert[fng] + 1;
	} // END (i; loop over neighbors of deleted particle)
	vert.clear();

	  // std::cout << "\n";     //////////////////////////////////////////////////////////////////////////////////////////////
	  /* for (i = 0; i < neigh_del.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << neigh_del[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sr_del.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << sr_del[2 * i] << " " << sr_del[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < sio_del.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sio_del[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sap_del.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sap_del[i] << " " << sap_del[i + 1] << " " << sap_del[i + 2] << " ; ";
	  } std::cout << " \n"; */

	con.put(del, nx, ny, nz);					// add new particle
												// co kdyz je nektera z nx, ny, nz mimo okno? put si asi diky periodicite poradi, ale poradi si i zbytek???
												//   konkretne poradi si fce face_dist? ne fci face_dist se musi predat pozice v oknë!!!!!!!!!!!
												// je potreba premapovat nx,ny,nz do okna ---> implementovano jiz ve fci bdma_step

	find_pos(ijk_add, q_add, del, &con);			// determine position of new particle
	con.compute_cell(c_add, ijk_add, q_add);	// compute cell for new added particle
	c_add.neighbors(neigh_add);					// compute its neighbors	
												//n_add = neigh_add.size();					// determine number of neighbors of added particle      !!!mozna zbytecne
	c_add.face_vertices(vert);					// {F,1PER} computes list of vertices for each face of the cell
	fng = 0;								// {F,1PER}

	// std::cout << ijk_add << " " << q_add << "\n"; ///////////////////////////////////////////////////////////////////////////

	k1_add = 0;
	
	vol = c_add.volume();
	part1 = 0;

	for (i = 0; i < neigh_add.size(); i++) {	// loop over the neighbors
		find_pos(ijk, q, neigh_add[i], &con);	// find position of the i-th neighbor

		// std::cout << con.p[ijk][3 * q] << " " << con.p[ijk][3 * q + 1] << " " << con.p[ijk][3 * q + 2] << " ; "; /////////////////

		face_dist(fng, vert, nx, ny, nz, xn, yn, zn, c_add);	// {F,PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
		xn = nx + 2 * xn; yn = ny + 2 * yn; zn = nz + 2 * zn;	// {F,PER} urcim skutecne souradnice tohoto souseda za stenou

		sr_add.push_back(ijk); sr_add.push_back(q); // add this neighbor to the second particle vector

		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {PER} nedoslo k preklopeni
			//part1 = part1 + V2(c_add, d);			// {1} add the value of V2 function to the first sum
			sio_add.push_back(0);					// neighbor is "in"
		}
		else {
			//part1 = part1 + V2(c_add, d, nx, ny, nz, xn, yn, zn);
			k1_add++;
			sio_add.push_back(k1_add);					// neighbor is k-th "out"
			sap_add.push_back(xn); sap_add.push_back(yn); sap_add.push_back(zn);  // storing real position of this neighbor
		}
		fng = fng + vert[fng] + 1;			// {1PER} set actual position in vector of face vertices
	} // END for (i; loop over neighbors of added particle)

	  // std::cout << "\n";    ////////////////////////////////////////////////////////////////////////////////////////////////////

	  /* std::cout << "add \n";     //////////////////////////////////////////////////////////////////////////////////////////////
	  for (i = 0; i < neigh_add.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << neigh_add[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sr_add.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << sr_add[2 * i] << " " << sr_add[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < sio_add.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sio_add[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sap_add.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sap_add[i] << " " << sap_add[i + 1] << " " << sap_add[i + 2] << " ; ";
	  } std::cout << " \n"; */

//	std::cout << "part1: " << part1 << " \n";  ////////////////////////////////////////////////////////////////////////////////

	// B - merge datovych struktur
	k = neigh_add.size();
	merge(neigh_add, neigh_del, sr_add, sr_del, k1_add, sio_add, sio_del, sap_add, sap_del);
	// slouceni se provede do vektoru sr_add, sio_add, sap_add
	// vektory sr_del, sio_del, sap_del se mi jeste budou hodit (jeste nebyla overena pripustnost pro vsechny zmenene castice)
	// vektor neigh_del byl podmnozinou vektoru neigh_add, pokud stale plati k == neigh_add.size() i po slouceni (tj. nebylo nic pridano)

	/* std::cout << "merge \n"; /////////////////////////////////////////////////////////////////////////////////////////////////
	for (i = 0; i < sr_add.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	std::cout << sr_add[2 * i] << " " << sr_add[2 * i + 1] << " ; ";
	} std::cout << " \n";
	for (i = 0; i < sio_add.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	std::cout << sio_add[i] << " ";
	} std::cout << " \n";
	for (i = 0; i < sap_add.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	std::cout << sap_add[i] << " " << sap_add[i + 1] << " " << sap_add[i + 2] << " ; ";
	} std::cout << " \n"; */

	// jiny pristup feasibility:
	cells = sr_add;	
	cells.push_back(ijk_add);
	cells.push_back(q_add);

//	std::cout << "NEW Feasibility: ";
	if (feasibility(con, cells, alfa, beta, B)) { /*std::cout << "OK. \n";*/ }
	else { erase(ijk_add, q_add, &con); con.put(del, x_del, y_del, z_del); return 0; }

	

	part1 = 0;
	part2 = 0;
	part3 = 0;
	l = 0;

	// 4.12.2017: domnenka, ze je nutne proverit feasibilitu i pro celou strukturu sr_add, a ne jenom pro sr_del - provede se spolecne s vypoctem part2
	// 24.10.2017: domnenka, ze i kdyz jsou neigh_add a neigh_del identicke, prece jen muze ve 3D dojit ke zmene sousedske struktury mezi nimi (???) 
	//				tedy feasibilita pro sr_del se musi overit vzdy a (ne)lze tak provest zaroven s vypoctem part2 a part3
	//if (k < neigh_add.size()) {   // nejsou-li neigh_add, neigh_del identicke musime overit feasibilitu i pro smazanou castici
	
	//std::cout << "TEST (move): k ... " << k << "\n";
	for (i = 0; i < k; i++) {				// {3} loop over the secondary particles
		ijk_i = sr_add[2 * i]; q_i = sr_add[2 * i + 1];
		con.compute_cell(d, ijk_i, q_i);  // {3} compute the cell of this primary neighbor

		if (sio_add[i] == 0) {
			x = con.p[ijk_i][3 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][3 * q_i + 1];
			z = con.p[ijk_i][3 * q_i + 2];
		}
		else {
			x = sap_add[3 * (sio_add[i] - 1)];
			y = sap_add[3 * (sio_add[i] - 1) + 1];
			z = sap_add[3 * (sio_add[i] - 1) + 2];
		}

		part1 = part1 + V2(c_add, d, nx, ny, nz, x, y, z);
		
		//std::cout << i << " " << sio_add[i] << " ... " << part1 << "\n";
	}

	k2 = 0;

	for (i = 0; i < sr_add.size() / 2; i++) {								// {2} loop over secondary particles 
		ijk_i = sr_add[2 * i]; q_i = sr_add[2 * i + 1];
		con.compute_cell(c, ijk_i, q_i);				// {2} compute cell of i-th secondary particle
		

		if (sio_add[i] == 0) {
			x = con.p[ijk_i][3 * q_i];					// computes coordinates of i-th particle
			y = con.p[ijk_i][3 * q_i + 1];
			z = con.p[ijk_i][3 * q_i + 2];
		}
		else {
			x = sap_add[3 * (sio_add[i] - 1)];
			y = sap_add[3 * (sio_add[i] - 1) + 1];
			z = sap_add[3 * (sio_add[i] - 1) + 2];
		}

		for (j = i + 1; j < sr_add.size() / 2; j++) {                   // {2} loop over secondary particles with "higher order" (to avoid double counting) 
			if (are_neighbors(c, sr_add[2 * j], sr_add[2 * j + 1], &con)) {	// {2} sousedi-li bunky i a j ... (i.e. if c_i ~ c_j)
				ijk_j = sr_add[2 * j]; q_j = sr_add[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		// {2} compute cell of the j-th secondary particle

				if (sio_add[j] == 0) {
					xn = con.p[ijk_j][3 * q_j];					// computes true coordinates of j-th particle
					yn = con.p[ijk_j][3 * q_j + 1];
					zn = con.p[ijk_j][3 * q_j + 2];
				}
				else {
					xn = sap_add[3 * (sio_add[j] - 1)];
					yn = sap_add[3 * (sio_add[j] - 1) + 1];
					zn = sap_add[3 * (sio_add[j] - 1) + 2];
				}

				part2 = part2 + V2(c, d, x, y, z, xn, yn, zn);
			
			} // END if (are_neighbors)
		} // END for (j; loop over secondary particles with higher order)
		c.neighbors(neigh);									           // {3} compute its neighbors ("secondary")
		c.face_vertices(vert);					// {3PER} computes list of vertices for each face of the cell
		fng = 0;								// {3PER}


		for (j = 0; j < neigh.size(); j++) {				           // {3} loop over neighbors of i-th secondary particle
			if (terciary(neigh[j], del, sr_add, &con)) {				 // verify if the neighbor is terciary or not (tj. nepatri mezi sekundarni castice a zaroven se take nejedna o nove pridanou castici
				find_pos(ijk, q, neigh[j], &con);			             // find its position
				tr.push_back(ijk); tr.push_back(q);			             // store this secondary neighbor to the vector of terciary particles 

																		 // pro volani fce face_dist musim pouzit stred bunky, tj souradnice ktere jsou v okne !!! a nedbat na to, zda jsou spravne
				face_dist(fng, vert, con.p[ijk_i][3 * q_i], con.p[ijk_i][3 * q_i + 1], con.p[ijk_i][3 * q_i + 2], xn, yn, zn, c);		// {3PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
																																																		// dale jiz ale mohu pouzit skutecne souradnice bodu bunky (x,y,z):
				xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;// {3PER} urcim skutecne souradnice tohoto souseda za stenou

																  //					if (xn == con.p[ijk][3 * q] && yn == con.p[ijk][3 * q + 1] && zn == con.p[ijk][3 * q + 2]) { // {3PER} nedoslo k preklopeni terciarni castice (terciarni castice je "in")
				if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {3PER} nedoslo k preklopeni
					tio.push_back(0);
					con.compute_cell(d, ijk, q);
					part3 = part3 + V2(c, d, x, y, z, xn, yn, zn);
					
				}
				else {								// terciarni castice je preklopena
					k2++; tio.push_back(k2);
					tap.push_back(xn); tap.push_back(yn); tap.push_back(zn); // storing real position of this neighbor
					con.compute_cell(d, ijk, q);
					part3 = part3 + V2(c, d, x, y, z, xn, yn, zn);
					
				} // END if..else (terciary particle is IN vs OUT)
				l++;
			} // END if (terciary)
			fng = fng + vert[fng] + 1;
		} // END for (j; loop over neigbors of i-th secondary particle)
		ntr.push_back(l);
	} // END for (i; loop over secondary particles)

//	std::cout << "part2: " << part2 << " \n";  /////////////////////////////////////////////////////////////////////
//	std::cout << "part3: " << part3 << " \n";  /////////////////////////////////////////////////////////////////////

											   /* for (i = 0; i < ntr.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
											   std::cout << ntr[i] << " ";
											   } std::cout << " \n";
											   for (i = 0; i < tr.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
											   std::cout << tr[2 * i] << " " << tr[2 * i + 1] << " ; ";
											   } std::cout << " \n";
											   for (i = 0; i < tio.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
											   std::cout << tio[i] << " ";
											   } std::cout << " \n";
											   for (i = 0; i < tap.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
											   std::cout << tap[i] << " " << tap[i + 1] << " " << tap[i + 2] << " ; ";
											   } std::cout << " \n"; */

	erase(ijk_add, q_add, &con);			// delete new particle
	con.put(del, x_del, y_del, z_del);		// add back old particle

	find_pos(ijk_del, q_del, del, &con);	// actualize position of old particle (ijk by se zmenit nemelo, ale q ano - optimalizovat!
	con.compute_cell(c_del, ijk_del, q_del);// compute the cell of deleted particle

	part4 = 0;

	for (i = 0; i < sr_del.size() / 2; i++) {					// {4} loop over secondary particles of DEL
		ijk_i = sr_del[2 * i]; q_i = sr_del[2 * i + 1];
		con.compute_cell(c, ijk_i, q_i);  // {4} compute the cell of this primary neighbor

		if (sio_del[i] == 0) {
			x = con.p[ijk_i][3 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][3 * q_i + 1];
			z = con.p[ijk_i][3 * q_i + 2];
		}
		else {
			x = sap_del[3 * (sio_del[i] - 1)];
			y = sap_del[3 * (sio_del[i] - 1) + 1];
			z = sap_del[3 * (sio_del[i] - 1) + 2];
		}

		part4 = part4 + V2(c_add, d, x_del, y_del, z_del, x, y, z);
		
	}
//	std::cout << "part4: " << part4 << " \n";  ////////////////////////////////////////////////////////////////////////

											   // std::cout << "k1_add, k1_del, k2: " << k1_add << " " << k1_del << " " << k2 << " \n"; //////////////////////////////////

	part5 = 0;
	part6 = 0;
	for (i = 0; i < sr_add.size() / 2; i++) {                              // {5,6} loop over secondary particles 
		ijk_i = sr_add[2 * i]; q_i = sr_add[2 * i + 1];
		if (sio_add[i] == 0) {
			x = con.p[ijk_i][3 * q_i];					// computes coordinates of i-th particle
			y = con.p[ijk_i][3 * q_i + 1];
			z = con.p[ijk_i][3 * q_i + 2];
		}
		else {
			x = sap_add[3 * (sio_add[i] - 1)];
			y = sap_add[3 * (sio_add[i] - 1) + 1];
			z = sap_add[3 * (sio_add[i] - 1) + 2];
		}
		con.compute_cell(c, ijk_i, q_i);			// {5,6} compute the cell of i-th primary neighbor

		for (j = i + 1; j < sr_add.size() / 2; j++) {				  // {5} loop over secondary particles with "higher order" (to avoid double counting)  	
			if (are_neighbors(c, sr_add[2 * j], sr_add[2 * j + 1], &con)) {	  // if i and j are neighbors... (i.e. if c_i ~ c_j)
				ijk_j = sr_add[2 * j]; q_j = sr_add[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		  // compute cell j

				if (sio_add[j] == 0) {
					xn = con.p[ijk_j][3 * q_j];					// computes coordinates of j-th particle
					yn = con.p[ijk_j][3 * q_j + 1];
					zn = con.p[ijk_j][3 * q_j + 2];
				}
				else {
					xn = sap_add[3 * (sio_add[j] - 1)];
					yn = sap_add[3 * (sio_add[j] - 1) + 1];
					zn = sap_add[3 * (sio_add[j] - 1) + 2];
				}

				part5 = part5 + V2(c, d, x, y, z, xn, yn, zn);
			}
		}

		if (i == 0) {
			for (j = 0; j < ntr[0]; j++) {							// {6} i=0; loop over its terciary particles		
				ijk_j = tr[2 * j]; q_j = tr[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		  // compute cell of this terciary particle

				if (tio[j] == 0) {
					xn = con.p[ijk_j][3 * q_j];								// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][3 * q_j + 1];
					zn = con.p[ijk_j][3 * q_j + 2];
				}
				else {
					xn = tap[3 * (tio[j] - 1)];
					yn = tap[3 * (tio[j] - 1) + 1];
					zn = tap[3 * (tio[j] - 1) + 2];
				}

				part6 = part6 + V2(c, d, x, y, z, xn, yn, zn);
			}
		}
		else {
			for (j = 0; j < ntr[i] - ntr[i - 1]; j++) {				// {6} i>0; loop over its terciary particles
				k = ntr[i - 1] + j;
				ijk_j = tr[2 * k]; q_j = tr[2 * k + 1];
				con.compute_cell(d, ijk_j, q_j);

				if (tio[k] == 0) {
					xn = con.p[ijk_j][3 * q_j];								// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][3 * q_j + 1];
					zn = con.p[ijk_j][3 * q_j + 2];
				}
				else {
					xn = tap[3 * (tio[k] - 1)];
					yn = tap[3 * (tio[k] - 1) + 1];
					zn = tap[3 * (tio[k] - 1) + 2];
				}

				part6 = part6 + V2(c, d, x, y, z, xn, yn, zn);
				
			}
		} // END if..else (i=0 vs i>0)

		
	} // END for (i; loop over secondary particles)

//	std::cout << "part5: " << part5 << " \n";  /////////////////////////////////////////////////////////////////////
//	std::cout << "part6: " << part6 << " \n";  /////////////////////////////////////////////////////////////////////

//	std::cout << "M: part1: " << -part1 << " | part2: " << -part2 << " | part3: " << -part3 << " | part4: " << part4 << " | part5: " << part5 << " | part6: " << part6 << " \n";

	return exp(theta*(-part1 - part2 - part3 + part4 + part5 + part6));
}


// ===============================================================================================================================================
// ===============================================================================================================================================
// ===============================================================================================================================================
// ===============================================================================================================================================
// ===============================================================================================================================================
// ===============================================================================================================================================
// ===============================================================================================================================================
// ===============================================================================================================================================
// ===============================================================================================================================================
// ===============================================================================================================================================
// ===============================================================================================================================================
// ===============================================================================================================================================

// ===============================================================================================================================================
// ==================================================================Laguerre=====================================================================

// 
double try_add(int id, double &nx, double &ny, double &nz, double &nrad, container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, const double &iota)
{
	// [in]		id;x,y,z; rad	the ID number, coordinates and radius of the particle.
	// [in]		con				the container with stored particles.
	// [in]		theta			parameter
	// [in]		alfa,beta,B		hardcore parameters.
	// [in]		iota			hardcore parameter of overlapping.

	unsigned int i, j, k, k1, k2, l, fng;
	int ijk, q, ijk_add, q_add, ijk_i, q_i, ijk_j, q_j;
	double part1, part2, part3, part4, part5;
	double xn, yn, zn, xnn, ynn, znn;
	int t1, t2, t3;

	std::vector<int> sr, tr;		// vectors storing ijk,q information about secondary, terciary particles
	std::vector<int> ntr;			// vector containing numbers of terciary particles for all secondary particles
	std::vector<int> sio, tio;		// information in/out for secondary, terciary particles
	std::vector<double> sap, tap;	// actual positions of secondary, terciary particles which are "out"
	std::vector<int> neigh, neighi;	// vector containing neighbor's ID
	std::vector<int> vert, verti;
	std::vector<int> cells;
	voronoicell_neighbor c, d;		// bunka s informaci o sousedech
	voronoicell di;

	//double alfa, beta; // asi je treba predat teto funkci parametricky
	//alfa = 0.02;
	//beta = 0.2;

	con.put(id, nx, ny, nz, nrad);   // pridani nove castice do containeru  
	find_pos(ijk_add, q_add, id, &con); // najde pozici nove castice 
										//
	// lepsi overit feasibility spolecne s vypoctem part1 (vypocet ma spolecne prvky)
	// Addel: pokus zda oddeleni techto castic nezvysi rychlost simulace (ve vetsine pripadu je pridani bodu zastaveno pripustnosti)

	// verifying the constrain on overlapping generators do not need computation of tessellation
	if (overlap_f(con, iota)) {}
	else { erase(ijk_add, q_add, &con); return 0; }
	// if iota >=0 then feasibility can be verified on the neighbourhood (cells) only, otherwise loose of the control threatens

	con.compute_cell(c, ijk_add, q_add);	// {1,F} compute cell for new added particle

	c.neighbors(neigh);						// {1,F} compute its neighbors			
											//c.face_vertices(vert);					// {1PER} computes list of vertices for each face of the cell
											//fng = 0;								// {1PER}


	k1 = 0;
	// spocti pravdepodobnost s kterou pridame castici; nejdrive spoctem vnitrek exponenciely po castech:
	// part1, part2, part3  se tykaji bunek po pridani nove castice
	//part1 = 0; // V2 pro nove pridanou bunku a jeji primarni sousedy (secondary particles)
	//fng = 0;

	// A  - ulozeni struktury sekundarnich castic
	for (i = 0; i < neigh.size(); i++) {	// {1,F} loop over the neighbors
		find_pos(ijk, q, neigh[i], &con);   // {1,F} find position of the i-th neighbor
											// std::cout << con.p[ijk][4 * q] << " " << con.p[ijk][4 * q + 1] << " " << con.p[ijk][4 * q + 2] << " ; "; /////////////////

		sr.push_back(ijk); sr.push_back(q); // add this neighbor to the second particle vector

		xn = con.p[ijk][4 * q]; yn = con.p[ijk][4 * q + 1]; zn = con.p[ijk][4 * q + 2];
		real_coo(nx, ny, nz, xn, yn, zn);

		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		 // {1PER} je-li pozice souseda v okne (tj. nedoslo k preklopeni)
																			 // std::cout << "IN \n"; ///////////////////////////////////////////////////////////////////////////////////
																			 //part1 = part1 + V2(c, d);    // {1} add the value of V2 function to the first sum
			sio.push_back(0);					// neighbor is "in"
		}
		else {
			//part1 = part1 + V2(c, d, false, con.p[ijk][4 * q], con.p[ijk][4 * q + 1], con.p[ijk][4 * q + 2], xn, yn, zn);
			k1++;
			sio.push_back(k1);					// neighbor is k-th "out"
			sap.push_back(xn); sap.push_back(yn); sap.push_back(zn);  // storing real position of this neighbor
		}
		// std::cout << V2(theta, c, d) << " \n"; ///////////////////////////////////////////////////////////////////////
		//fng = fng + vert[fng] + 1;			// {1PER} set actual position in vector of face vertices
	}

	// B - overeni feasibility
	cells = sr;
	cells.push_back(ijk_add); cells.push_back(q_add);

	//	1) Feasibility on the neighborhood pr+sr:
		if (feasibility(con, cells, alfa, beta, B)) { /*std::cout << "OK. \n";*/ }
		else { erase(ijk_add, q_add, &con); return 0; }
	// PROBLEM: overit feasibility pro cells pouze je nedostatecne, chtelo by to nalezt minimalni mnozinu 
	// solution: pokud iota >= 0 pak je dostatecne

	//	2) Feasibility of the whole container:
	//	if (feasibility(con, alfa, beta, B)) { /*std::cout << "OK. \n";*/ } else { erase(ijk_add, q_add, &con); return 0; }

	//	3) Feasibility on the subwindow:
	//if (feasibility_subwindows(con, alfa, beta, B, nx, ny, nz, nrad)) { /*std::cout << "OK. \n";*/ }
	//else { erase(ijk_add, q_add, &con); return 0; }


	part1 = 0;
	for (i = 0; i < sr.size() / 2; i++) {
		ijk_i = sr[2 * i]; q_i = sr[2 * i + 1];
		con.compute_cell(d, ijk_i, q_i);

		if (sio[i] == 0) {
			xn = con.p[ijk_i][4 * q_i];					// {F} computes coordinates of i-th particle
			yn = con.p[ijk_i][4 * q_i + 1];
			zn = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			xn = sap[3 * (sio[i] - 1)];
			yn = sap[3 * (sio[i] - 1) + 1];
			zn = sap[3 * (sio[i] - 1) + 2];
		}

		part1 = part1 + V2(c, d, nx, ny, nz, xn, yn, zn);
		/*
		if (sio[i] == 0) {
		part1 = part1 + V2(c, d);
		}
		else {
		part1 = part1 + V2(c, d, false, con.p[ijk_i][4 * q_i], con.p[ijk_i][4 * q_i + 1], con.p[ijk_i][4 * q_i + 2],
		sap[3 * (sio[i] - 1)], sap[3 * (sio[i] - 1) + 1], sap[3 * (sio[i] - 1) + 2]);
		}*/
	}



	/*std::cout << " \n"; /////////////////////////////////////////////////////////////////////////////////////////////////
	for (i = 0; i < sr.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	std::cout << sr[2 * i] << " " << sr[2 * i + 1] << " ; ";
	} std::cout << " \n";
	for (i = 0; i < sio.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	std::cout << sio[i] << " ";
	} std::cout << " \n";
	for (i = 0; i < sap.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	std::cout << sap[i] << " " << sap[i + 1] << " " << sap[i + 2] << " ; ";
	} std::cout << " \n"; */

	//	std::cout << "part1: " << part1 << " \n"; ///////////////////////////////////////////////////////////////////////////
	// std::cout << sr.size() << " sr \n";  /////////////////////////////////////////////////////////////////////////////
	// for (i = 0; i < sr.size() / 2; i++) { std::cout << sr[2 * i] << " " << sr[2 * i + 1] << " \n"; } //////////////////

	k2 = 0; l = 0;
	part2 = 0; // V2 pro primarni sousedy (po pridani)
	part3 = 0; // V2 pro primarni a sekundarni sousedy (po pridani)
	for (i = 0; i < sr.size() / 2; i++) {                           // {2,3} loop over secondary particles
		ijk_i = sr[2 * i]; q_i = sr[2 * i + 1];

		if (sio[i] == 0) {
			xnn = con.p[ijk_i][4 * q_i];							// {2,3PER} computes coordinates of i-th particle
			ynn = con.p[ijk_i][4 * q_i + 1];
			znn = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			xnn = sap[3 * (sio[i] - 1)];
			ynn = sap[3 * (sio[i] - 1) + 1];
			znn = sap[3 * (sio[i] - 1) + 2];
		}

		con.compute_cell(c, ijk_i, q_i);				// {2,3} compute the cell of the i-th secundary particle

		for (j = i + 1; j < sr.size() / 2; j++) {					// {2} loop over secondary particles with "higher order" (to avoid double counting)  	
			ijk_j = sr[2 * j]; q_j = sr[2 * j + 1];
			if (are_neighbors(c, ijk_j, q_j, &con)) {	  // if i and j are neighbors... (i.e. if c_i ~ c_j)
				con.compute_cell(d, ijk_j, q_j);		  // compute cell j

				if (sio[j] == 0) {
					xn = con.p[ijk_j][4 * q_j];							// {2,3PER} computes coordinates of i-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = sap[3 * (sio[j] - 1)];
					yn = sap[3 * (sio[j] - 1) + 1];
					zn = sap[3 * (sio[j] - 1) + 2];
				}

				part2 = part2 + V2(c, d, xnn, ynn, znn, xn, yn, zn);
			}
		}

		c.neighbors(neigh);											// {3} compute its neighbors
																	//c.face_vertices(vert);										// {3PER} computes list of vertices for each face of the cell
																	//fng = 0;													// {3PER}

		for (j = 0; j < neigh.size(); j++) {						// {3} loop over these neighbors

			if (terciary(neigh[j], id, sr, &con)) {			          // {3} je-li soused "terciarni castice", tj je sousedem souseda pridane castice a zaroven neni sousedem pridane castice
				find_pos(ijk, q, neigh[j], &con);					    // najdi pozici souseda
				tr.push_back(ijk); tr.push_back(q);					    // uloz souseda mezi terciarni 

				real_coo(xnn, ynn, znn, xn, yn, zn);

				if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {  // {3PER} nedoslo k preklopeni (terciarni castice je "in")
					tio.push_back(0);
					con.compute_cell(d, ijk, q);
					part3 = part3 + V2(c, d, xnn, ynn, znn, xn, yn, zn);
				}
				else {								// doslo k preklopeni (terciarni catice je "out")
					k2++; tio.push_back(k2);
					tap.push_back(xn); tap.push_back(yn); tap.push_back(zn); // storing real position of this neighbor

					con.compute_cell(d, ijk, q);
					part3 = part3 + V2(c, d, xnn, ynn, znn, xn, yn, zn);
					/*part3 = part3 + V2(c, d, false, con.p[ijk][3 * q], con.p[ijk][3 * q + 1], con.p[ijk][3 * q + 2],
					xn, yn, zn);*/
					// con.p[sr[2 * i]][3 * sr[2 * i + 1]], con.p[sr[2 * i]][3 * sr[2 * i + 1] + 1], con.p[sr[2 * i]][3 * sr[2 * i + 1] + 2], xn, yn, zn);

				} // END if..else (preklopeni)
				l++;
			} // END if (terciary)
			  //fng = fng + vert[fng] + 1;			// {2,3PER} set actual position in vector of face vertices (posun na dalsi stenu)
		} // END for (j; cyklus over neighbors of secondary particle i)
		ntr.push_back(l);									// uloz pocet terciarnich castic souseda i
	} // END for (i; cyklus over sr)
	  //	std::cout << "part2: " << part2 << " \n";  ///////////////////////////////////////////////////////////////////// 


	  // std::cout << terciary(5, id, sr, &con) << " \n"; ///////////////////////////////////////////////////////////////
	  //	std::cout << "part3: " << part3 << " \n";  /////////////////////////////////////////////////////////////////////

	  /*std::cout << tr.size() << " tr \n";  ///////////////////////////////////////////////////////////////////////////
	  std::cout << ntr.size() << " ntr \n"; //////////////////////////////////////////////////////////////////////////
	  for (i = 0; i < ntr.size(); i++) {
	  std::cout << ntr[i] << " ";
	  } std::cout << " \n";  ///////////////////////////////////////////////////////////////////////////////////////////////
	  for (i = 0; i < tr.size() / 2; i++) {
	  std::cout << tr[2 * i] << " " << tr[2 * i + 1] << " ; ";
	  } std::cout << " \n";  /////////////////////////////////////////////////////////////////////////////////////////
	  for (i = 0; i < tio.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << tio[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < tap.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << tap[i] << " " << tap[i + 1] << " " << tap[i + 2] << " ; ";
	  } std::cout << " \n"; */
	  // std::cout << con.total_particles() << " \n"; ///////////////////////////////////////////////////////////////////

	  // nyni jsme spocetli cast vnitrku exponenciely tykajici se stavu po pridani bunky;
	  // dale zbyva cast tykajici se stavu pred pridanim 
	  // - odstranim novou castici a vyuziju jiz spoctenych sekundarnich a terciarnich castic
	  // - nevyhodou je ze musim znovu cyklit, toho by se slo zbavit tim, ze bych mel ulozene dva containery - pro stavy
	  //	 pred a po  (je dulezitejsi casova uspora oproti pametove?)

	erase(ijk_add, q_add, &con); // smaz novou castici at muzeme spocitat stav pred pridanim

								 // std::cout << con.total_particles() << " \n";  //////////////////////////////////////////////////////////////////

								 // part4, part5  se tykaji bunek pred pridanim nove castice
	part4 = 0; // V2 pro primarni sousedy (pred pridanim)  // cast je totozna s part2
	part5 = 0; // V2 pro primarni a sekundarni sousedy (pred pridanim)
	for (i = 0; i < sr.size() / 2; i++) {							// {4,5} loop over secondary particles
		ijk_i = sr[2 * i]; q_i = sr[2 * i + 1];
		con.compute_cell(c, ijk_i, q_i);							// {4,5} compute its cells
																	// c.neighbors(neigh);										// {4,5} compute its neighbors
		if (sio[i] == 0) {
			xn = con.p[ijk_i][4 * q_i];								// {F} computes coordinates of i-th particle
			yn = con.p[ijk_i][4 * q_i + 1];
			zn = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			xn = sap[3 * (sio[i] - 1)];
			yn = sap[3 * (sio[i] - 1) + 1];
			zn = sap[3 * (sio[i] - 1) + 2];
		}

		for (j = i + 1; j < sr.size() / 2; j++) {					// {4} loop over secondary particles with "higher order" (to avoid double counting)  	
			ijk_j = sr[2 * j]; q_j = sr[2 * j + 1];
			if (are_neighbors(c, ijk_j, q_j, &con)) {	  // if i and j are neighbors... (i.e. if c_i ~ c_j)
				con.compute_cell(d, ijk_j, q_j);

				if (sio[j] == 0) {
					xnn = con.p[ijk_j][4 * q_j];								// {F} computes coordinates of i-th particle
					ynn = con.p[ijk_j][4 * q_j + 1];
					znn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xnn = sap[3 * (sio[j] - 1)];
					ynn = sap[3 * (sio[j] - 1) + 1];
					znn = sap[3 * (sio[j] - 1) + 2];
				}

				part4 = part4 + V2(c, d, xn, yn, zn, xnn, ynn, znn);
			}
		}

		if (i == 0) {
			for (j = 0; j < ntr[0]; j++) {							// {5} i=0; loop over its terciary particles
				ijk_j = tr[2 * j]; q_j = tr[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		  // compute cell of this terciary particle

				if (tio[j] == 0) {
					xnn = con.p[ijk_j][4 * q_j];								// {F} computes coordinates of i-th particle
					ynn = con.p[ijk_j][4 * q_j + 1];
					znn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xnn = tap[3 * (tio[j] - 1)];
					ynn = tap[3 * (tio[j] - 1) + 1];
					znn = tap[3 * (tio[j] - 1) + 2];
				}

				part5 = part5 + V2(c, d, xn, yn, zn, xnn, ynn, znn);
			}
		}
		else {
			for (j = 0; j < ntr[i] - ntr[i - 1]; j++) {				// {5} i>0; loop over its terciary particles
				k = ntr[i - 1] + j;
				ijk_j = tr[2 * k]; q_j = tr[2 * k + 1];
				con.compute_cell(d, ijk_j, q_j);

				if (tio[k] == 0) {
					xnn = con.p[ijk_j][4 * q_j];								// {F} computes coordinates of i-th particle
					ynn = con.p[ijk_j][4 * q_j + 1];
					znn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xnn = tap[3 * (tio[k] - 1)];
					ynn = tap[3 * (tio[k] - 1) + 1];
					znn = tap[3 * (tio[k] - 1) + 2];
				}

				part5 = part5 + V2(c, d, xn, yn, zn, xnn, ynn, znn);



			}
		} // END if..else (i==0 vs i>0)

	} // END for (i; cyklus over secondary particles)

	  //	std::cout << "part4: " << part4 << " \n";  ///////////////////////////////////////////////////////////////////// 
	  //	std::cout << "part5: " << part5 << " \n";  ///////////////////////////////////////////////////////////////////// 

	  //	std::cout << "A: part1: " << -part1 << " | part2: " << -part2 << " | part3: " << -part3 << " | part4: " << part4 << " | part5: " << part5 << " \n";

	return exp(theta*(-part1 - part2 - part3 + part4 + part5));
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
// _____________________________________________________________________________________________________________________________________________


// ---------------------------------------------------------------------------------------------------------------------------------------------
// _____________________________________________________________________________________________________________________________________________
// _____________________________________________________________________________________________________________________________________________



// 395 - 755

double try_delete(int &ijk_del, int &q_del, container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, const double &iota)
{
	// [in]		ijk_del,q_del	position of the particle.
	// [in]		con				the container with stored particles.
	// [in]		theta			parameter
	// [in]		alfa,beta,B		hardcore parameters.

	unsigned int i, j, k, k1, k2, l, fng;
	int ijk, q, ijk_i, q_i, ijk_j, q_j, del;
	double part1, part2, part3, part4, part5;
	double x_del, y_del, z_del, x, y, z, xn, yn, zn;
	double rad_del;

	std::vector<int> pr, sr, tr; // misto particle_order class snadnejsi pouzivat vektory pro ukladani ijk,q castic
	std::vector<int> sio, tio;
	std::vector<double> sap, tap;
	std::vector<int> ntr, cells;
	std::vector<int> neigh;  // vektor pro prirazeni ID sousedu
	std::vector<int> vert;
	voronoicell_neighbor c_del, c, d;  // bunka s informaci o sousedech

									   //double alfa, beta; // asi je treba predat teto funkci parametricky
									   //alfa = 0.05;
									   //beta = 0.5;

									   // priprava: najdu souradnice a ID mazane castice, spoctu jeji bunku, pak jeji sousedy a vrcholy sten
									   //		loop pres sousedy ulozi do sr pozice ijk,q techto sekundarnich castic, face_dist spocte skutecne polohy 
									   //		techto sousedu, do sio se ulozi informace o preklopeni a do sap pripadne skutecna poloha sekundarni castice

	x_del = con.p[ijk_del][4 * q_del];      // coordinates of deleted particle
	y_del = con.p[ijk_del][4 * q_del + 1];
	z_del = con.p[ijk_del][4 * q_del + 2];
	rad_del = con.p[ijk_del][4 * q_del + 3];
	del = con.id[ijk_del][q_del];			// ID of deleted particle

	con.compute_cell(c_del, ijk_del, q_del);// {} compute the cell of deleted particle
	c_del.neighbors(neigh);					// {} compute its neighbors (in following "primary")
											//c_del.face_vertices(vert);				// {PER} computes list of vertices for each face of the cell
											//fng = 0;								// {PER}

											//std::cout << del << " " << ijk_del << "\n";  /////////////////////////////////////////////////////////////////////////
											//std::cout << "uvod /n";  /////////////////////////////////////////////////////////////////////////////////////////////
											// nejdrive se spocetli udaje o mazane castici, jeji bunka, sousede, vrcholy sten; pak se castice smaze, a az pak se 
											//    ulozi vektory sr,sio,sap (udaje o castici c_del zustaly ulozeny, takze nevadi ze tato bunka uz nebude existovat)
											//    (nutne je toto poradi, protoze fce erase meni hodnoty q, tj pozice v boxu)

											// smazani castice a vypocet stavu po smazani (mezi prvnimi nas zajima feasibilita, bez ni nema smysl pokracovat)
	erase(ijk_del, q_del, &con); // smaz castici at muzeme spocitat stav po smazani, smazani zmeni strukturu

	// overlapping verification
	if (overlap_f(con, iota)) {}
	else { con.put(del, x_del, y_del, z_del, rad_del); return 0; }

								 //std::cout << "smazano \n"; ///////////////////////////////////////////////////////////////////////////////////////////

	k1 = 0;
	for (i = 0; i < neigh.size(); i++) {	// {} loop over the neighbors of deleted particle
		find_pos(ijk, q, neigh[i], &con);   // {} find position of the i-th primary neighbor after deleting 

		sr.push_back(ijk); sr.push_back(q); // add this primary neighbor to the vector of secondary particles

		xn = con.p[ijk][4 * q]; yn = con.p[ijk][4 * q + 1]; zn = con.p[ijk][4 * q + 2];
		real_coo(x_del, y_del, z_del, xn, yn, zn);

		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {PER} nedoslo k preklopeni
			sio.push_back(0);					// neighbor is "in"
		}
		else {
			k1++;
			sio.push_back(k1);					// neighbor is k-th "out"
			sap.push_back(xn); sap.push_back(yn); sap.push_back(zn);  // storing real position of this neighbor
		}

		//fng = fng + vert[fng] + 1;
	} // END for (i; cyklus over neighbors of deleted particle)

	  /* for (i = 0; i < neigh.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << neigh[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sr.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << sr[2 * i] << " " << sr[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < sio.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sio[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sap.size()/3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sap[i] << " " << sap[i + 1] << " " << sap[i + 2] << " ; ";
	  } std::cout << " \n"; */
	  //std::cout << "priprava /n"; //////////////////////////////////////////////////////////////////////////////////////////

	  // Addel: striktni oddeleni overeni pripustnosti od ostatnich ukonu, pripustnost je prvorada


	  // jiny pristup feasibility:
	cells = sr;

	//	1) Feasibility of the neighborhood pr+sr:
		if (feasibility(con, cells, alfa, beta, B)) { /*std::cout << "OK. \n";*/ }
		else { con.put(del, x_del, y_del, z_del, rad_del); return 0; }

	//	2) Feasibility of the whole container:
	//	if (feasibility(con, alfa, beta, B)) { /*std::cout << "OK. \n";*/ }
	//	else { con.put(del, x_del, y_del, z_del, rad_del); return 0; }

	//	3) Feasibility on the subwindow:
	//if (feasibility_subwindows(con, alfa, beta, B, x_del, y_del, z_del, rad_del)) { /*std::cout << "OK. \n";*/ }
	//else { con.put(del, x_del, y_del, z_del, rad_del); return 0; }

	//-----------------------------------------------------------------------------------



	// B  - vypocet struktury terciarnich castic a part 1, part 2 
	part1 = 0;
	part2 = 0;
	l = 0;
	k2 = 0;
	for (i = 0; i < sr.size() / 2; i++) {                           // {1,2,F} loop over secondary particles
		ijk_i = sr[2 * i]; q_i = sr[2 * i + 1];
		con.compute_cell(c, ijk_i, q_i);							// {1,2,F} compute cell of i-th secondary particle

		if (sio[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap[3 * (sio[i] - 1)];
			y = sap[3 * (sio[i] - 1) + 1];
			z = sap[3 * (sio[i] - 1) + 2];
		}

		for (j = i + 1; j < sr.size() / 2; j++) {                   // {1,F} loop over secondary particles with "higher order" (to avoid double counting) 
			if (are_neighbors(c, sr[2 * j], sr[2 * j + 1], &con)) {	// {1,F} sousedi-li bunky i a j ... (i.e. if c_i ~ c_j)
				ijk_j = sr[2 * j]; q_j = sr[2 * j + 1];
				// std::cout << "are \n"; ///////////////////////////////////////////////////////////////////////////////
				if (sio[j] == 0) {
					xn = con.p[ijk_j][4 * q_j];				// {F} computes coordinates of j-th neighbor
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = sap[3 * (sio[j] - 1)];
					yn = sap[3 * (sio[j] - 1) + 1];
					zn = sap[3 * (sio[j] - 1) + 2];
				}

				con.compute_cell(d, ijk_j, q_j);					// {1} compute cell of the j-th secondary particle

				part1 = part1 + V2(c, d, x, y, z, xn, yn, zn);

			} // END if (are_neighbors)
		} // END for (j; secondary particles with higher order - j>i)
		c.neighbors(neigh);						// {2} compute its neighbors ("secondary")
												//c.face_vertices(vert);					// {2PER} computes list of vertices for each face of the cell
												//fng = 0;								// {2PER}


		for (j = 0; j < neigh.size(); j++) {				           // {2} loop over terciary particles (secondary neighbors)
			if (terciary(neigh[j], sr, &con)) {							 // je-li castice terciarni, tj soused sekundarni castice, ktery zaroven neni sekundarni castici
				find_pos(ijk, q, neigh[j], &con);			             // find its position

				tr.push_back(ijk); tr.push_back(q);			             // store this secondary neighbor to the vector of terciary particles 

				real_coo(x, y, z, xn, yn, zn);
				//					if (xn == con.p[ijk][3 * q] && yn == con.p[ijk][3 * q + 1] && zn == con.p[ijk][3 * q + 2]) { // {2PER} nedoslo k preklopeni terciarni castice (terciarni castice je "in")
				if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {2PER} nedoslo k preklopeni terciarni castice (terciarni castice je "in")
					tio.push_back(0);
					con.compute_cell(d, ijk, q);
					part2 = part2 + V2(c, d, x, y, z, xn, yn, zn);


				}
				else {								// terciarni castice je preklopena
					k2++; tio.push_back(k2);
					tap.push_back(xn); tap.push_back(yn); tap.push_back(zn); // storing real position of this neighbor

					con.compute_cell(d, ijk, q);
					part2 = part2 + V2(c, d, x, y, z, xn, yn, zn);

				} // END if..else (terciary particle is IN vs OUT)
				l++;
			} // END if (terciary)
			  //fng = fng + vert[fng] + 1;
		}  // END for (j; loop over neighbors of i-th secondary particle)
		ntr.push_back(l);									         // {2} uloz pocet sekundarnich sousedu souseda i
	} // END for (i; loop over secondary particles)


	  //	std::cout << "part1: " << part1 << " \n";  //////////////////////////////////////////////////////////////////////
	  //	std::cout << "part2: " << part2 << " \n";  //////////////////////////////////////////////////////////////////////

	  /* for (i = 0; i < ntr.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << ntr[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < tr.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << tr[2 * i] << " " << tr[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < tio.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << tio[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < tap.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << tap[i] << " " << tap[i + 1] << " " << tap[i + 2] << " ; ";
	  } std::cout << " \n"; */


	con.put(del, x_del, y_del, z_del, rad_del);
	// castice se prida sice asi do stejneho boxu, ale rozhodne ne na sve puvodni misto (nybrz na konec) !!!!!!!!
	//    !!!!!!! je proto potreba urcit opravdovou polohu teto castice !!!!!!!!!!!!!
	find_pos(ijk_del, q_del, del, &con);

	part3 = 0;
	for (i = 0; i < sr.size() / 2; i++) {				// {3} loop over the secondary particles
		ijk_i = sr[2 * i]; q_i = sr[2 * i + 1];
		con.compute_cell(d, ijk_i, q_i);  // {3} compute the cell of this primary neighbor

		if (sio[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap[3 * (sio[i] - 1)];
			y = sap[3 * (sio[i] - 1) + 1];
			z = sap[3 * (sio[i] - 1) + 2];
		}

		part3 = part3 + V2(c_del, d, x_del, y_del, z_del, x, y, z);


	}
	//	std::cout << "part3: " << part3 << " \n";  /////////////////////////////////////////////////////////////////////

	part4 = 0;
	part5 = 0;
	for (i = 0; i < sr.size() / 2; i++) {                              // {4,5} loop over secondary particles    
		ijk_i = sr[2 * i]; q_i = sr[2 * i + 1];
		if (sio[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap[3 * (sio[i] - 1)];
			y = sap[3 * (sio[i] - 1) + 1];
			z = sap[3 * (sio[i] - 1) + 2];
		}
		con.compute_cell(c, ijk_i, q_i);				// {4,5} compute the cell of i-th primary neighbor   

		for (j = i + 1; j < sr.size() / 2; j++) {					// {4} loop over secondary particles with "higher order" (to avoid double counting)  	
			if (are_neighbors(c, sr[2 * j], sr[2 * j + 1], &con)) {	  // if i and j are neighbors... (i.e. if c_i ~ c_j)
				ijk_j = sr[2 * j]; q_j = sr[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);					  // compute cell j

				if (sio[j] == 0) {
					xn = con.p[ijk_j][4 * q_j];					// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = sap[3 * (sio[j] - 1)];
					yn = sap[3 * (sio[j] - 1) + 1];
					zn = sap[3 * (sio[j] - 1) + 2];
				}

				part4 = part4 + V2(c, d, x, y, z, xn, yn, zn);
			}
		}

		if (i == 0) {
			for (j = 0; j < ntr[0]; j++) {							// {5} i=0; loop over its terciary particles
				ijk_j = tr[2 * j]; q_j = tr[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		  // compute cell of this terciary particle
				if (tio[j] == 0) {
					xn = con.p[ijk_j][4 * q_j];								// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = tap[3 * (tio[j] - 1)];
					yn = tap[3 * (tio[j] - 1) + 1];
					zn = tap[3 * (tio[j] - 1) + 2];
				}

				part5 = part5 + V2(c, d, x, y, z, xn, yn, zn);
			}
		}
		else {
			for (j = 0; j < ntr[i] - ntr[i - 1]; j++) {				// {5} i>0; loop over its terciary particles
				k = ntr[i - 1] + j;
				ijk_j = tr[2 * k]; q_j = tr[2 * k + 1];
				con.compute_cell(d, ijk_j, q_j);

				if (tio[k] == 0) {
					xn = con.p[ijk_j][4 * q_j];								// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = tap[3 * (tio[k] - 1)];
					yn = tap[3 * (tio[k] - 1) + 1];
					zn = tap[3 * (tio[k] - 1) + 2];
				}

				part5 = part5 + V2(c, d, x, y, z, xn, yn, zn);
			}
		} // END if..else (i==0 vs i>0)

	}  // END for (i; loop over secondary particles)

	   //	std::cout << "part4: " << part4 << " \n";  /////////////////////////////////////////////////////////////////////
	   //	std::cout << "part5: " << part5 << " \n";  /////////////////////////////////////////////////////////////////////

	   //	std::cout << "D: part1: " << -part1 << " | part2: " << -part2 << " | part3: " << part3 << " | part4: " << part4 << " | part5: " << part5 << " \n";

	return exp(theta*(-part1 - part2 + part3 + part4 + part5));
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
// _____________________________________________________________________________________________________________________________________________


// ---------------------------------------------------------------------------------------------------------------------------------------------
// _____________________________________________________________________________________________________________________________________________
// _____________________________________________________________________________________________________________________________________________


// "MOVE" consists from move and change of radius; these two actions can be treated together or separatedly

double try_MOVE(int &ijk_del, int &q_del, double &nx, double &ny, double &nz, double &nrad, container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, const double &iota)
{
	// [in]		ijk_del,q_del	position of the particle to be deleted.
	// [in]		x,y,z,rad		coordinates and radius of the move proposal.  -- misto tohoto by stacilo predat hodnotu sigma
	// [in]		con				the container with stored particles.
	// [in]		theta			parameter
	// [in]		alfa,beta,B		hardcore parameters.

	unsigned int i, j, k, l, k1_add, k1_del, k2, fng;
	int ijk, q, ijk_add, q_add, del, ijk_i, q_i, ijk_j, q_j;
	double part1, part2, part3, part4, part5, part6;
	double x_del, y_del, z_del, x, y, z, xn, yn, zn;
	double vol, rad_del;

	std::vector<int> sr_add, sr_del, tr; // misto particle_order class snadnejsi pouzivat vektory pro ukladani ijk,q castic
	std::vector<int> sio_add, sio_del, tio;
	std::vector<double> sap_add, sap_del, tap;
	std::vector<int> ntr, cells;
	std::vector<int> neigh_del, neigh_add, neigh, neighi;  // vektor pro prirazeni ID sousedu
	std::vector<int> vert, verti;
	voronoicell_neighbor c_del, c_add, c, d;  // bunka s informaci o sousedech


	x_del = con.p[ijk_del][4 * q_del];      // coordinates of deleted particle
	y_del = con.p[ijk_del][4 * q_del + 1];
	z_del = con.p[ijk_del][4 * q_del + 2];
	rad_del = con.p[ijk_del][4 * q_del + 3];

	del = con.id[ijk_del][q_del];			// ID of deleted particle

	con.compute_cell(c_del, ijk_del, q_del);	// compute the cell of deleted particle
	c_del.neighbors(neigh_del);					// compute its neighbors (in following "primary")
												//n_del = neigh_del.size();					// determine number of neighbors of deleted particle
												//c_del.face_vertices(vert);					// {PER} computes list of vertices for each face of the cell
												//fng = 0;

	erase(ijk_del, q_del, &con);				// delete old particle
												// std::cout << "del \n";     //////////////////////////////////////////////////////////////////////////////////////////////

	k1_del = 0;
	for (i = 0; i < neigh_del.size(); i++) {
		find_pos(ijk, q, neigh_del[i], &con);

		sr_del.push_back(ijk); sr_del.push_back(q);
		// std::cout << con.p[ijk][4 * q] << " " << con.p[ijk][4 * q + 1] << " " << con.p[ijk][4 * q + 2] << " ; "; /////////////////

		xn = con.p[ijk][4 * q]; yn = con.p[ijk][4 * q + 1]; zn = con.p[ijk][4 * q + 2];
		real_coo(x_del, y_del, z_del, xn, yn, zn);

		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {PER} nedoslo k preklopeni
			sio_del.push_back(0);					// neighbor is "in"
		}
		else {
			k1_del++;
			sio_del.push_back(k1_del);					// neighbor is k-th "out"
			sap_del.push_back(xn); sap_del.push_back(yn); sap_del.push_back(zn);  // storing real position of this neighbor
		}

		//fng = fng + vert[fng] + 1;
	} // END (i; loop over neighbors of deleted particle)
	  //vert.clear();

	  // std::cout << "\n";     //////////////////////////////////////////////////////////////////////////////////////////////
	  /* for (i = 0; i < neigh_del.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << neigh_del[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sr_del.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << sr_del[2 * i] << " " << sr_del[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < sio_del.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sio_del[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sap_del.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sap_del[i] << " " << sap_del[i + 1] << " " << sap_del[i + 2] << " ; ";
	  } std::cout << " \n"; */

	con.put(del, nx, ny, nz, nrad);					// add new particle
														// co kdyz je nektera z nx, ny, nz mimo okno? put si asi diky periodicite poradi, ale poradi si i zbytek???
														//   konkretne poradi si fce face_dist? ne fci face_dist se musi predat pozice v oknë!!!!!!!!!!!
														// je potreba premapovat nx,ny,nz do okna ---> implementovano jiz ve fci bdma_step


	find_pos(ijk_add, q_add, del, &con);			// determine position of new particle

	if (overlap_f(con, iota)) {}
	else { erase(ijk_add, q_add, &con); con.put(del, x_del, y_del, z_del, rad_del); return 0; }


	con.compute_cell(c_add, ijk_add, q_add);	// compute cell for new added particle
	c_add.neighbors(neigh_add);					// compute its neighbors	
												//n_add = neigh_add.size();					// determine number of neighbors of added particle      !!!mozna zbytecne
												//c_add.face_vertices(vert);					// {F,1PER} computes list of vertices for each face of the cell
												//fng = 0;								// {F,1PER}

												// std::cout << ijk_add << " " << q_add << "\n"; ///////////////////////////////////////////////////////////////////////////

	k1_add = 0;

	vol = c_add.volume();
	part1 = 0;

	for (i = 0; i < neigh_add.size(); i++) {	// loop over the neighbors
		find_pos(ijk, q, neigh_add[i], &con);	// find position of the i-th neighbor

												// std::cout << con.p[ijk][4 * q] << " " << con.p[ijk][4 * q + 1] << " " << con.p[ijk][4 * q + 2] << " ; "; /////////////////

		sr_add.push_back(ijk); sr_add.push_back(q); // add this neighbor to the second particle vector

		xn = con.p[ijk][4 * q]; yn = con.p[ijk][4 * q + 1]; zn = con.p[ijk][4 * q + 2];
		real_coo(nx, ny, nz, xn, yn, zn);

		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {PER} nedoslo k preklopeni
																			//part1 = part1 + V2(c_add, d);			// {1} add the value of V2 function to the first sum
			sio_add.push_back(0);					// neighbor is "in"
		}
		else {
			//part1 = part1 + V2(c_add, d, nx, ny, nz, xn, yn, zn);
			k1_add++;
			sio_add.push_back(k1_add);					// neighbor is k-th "out"
			sap_add.push_back(xn); sap_add.push_back(yn); sap_add.push_back(zn);  // storing real position of this neighbor
		}
		//fng = fng + vert[fng] + 1;			// {1PER} set actual position in vector of face vertices
	} // END for (i; loop over neighbors of added particle)

	  // std::cout << "\n";    ////////////////////////////////////////////////////////////////////////////////////////////////////

	  /* std::cout << "add \n";     //////////////////////////////////////////////////////////////////////////////////////////////
	  for (i = 0; i < neigh_add.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << neigh_add[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sr_add.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << sr_add[2 * i] << " " << sr_add[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < sio_add.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sio_add[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sap_add.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sap_add[i] << " " << sap_add[i + 1] << " " << sap_add[i + 2] << " ; ";
	  } std::cout << " \n"; */

	  //	std::cout << "part1: " << part1 << " \n";  ////////////////////////////////////////////////////////////////////////////////

	  // B - merge datovych struktur
	k = neigh_add.size();
	merge(neigh_add, neigh_del, sr_add, sr_del, k1_add, sio_add, sio_del, sap_add, sap_del);
	// slouceni se provede do vektoru sr_add, sio_add, sap_add
	// vektory sr_del, sio_del, sap_del se mi jeste budou hodit (jeste nebyla overena pripustnost pro vsechny zmenene castice)
	// vektor neigh_del byl podmnozinou vektoru neigh_add, pokud stale plati k == neigh_add.size() i po slouceni (tj. nebylo nic pridano)

	/* std::cout << "merge \n"; /////////////////////////////////////////////////////////////////////////////////////////////////
	for (i = 0; i < sr_add.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	std::cout << sr_add[2 * i] << " " << sr_add[2 * i + 1] << " ; ";
	} std::cout << " \n";
	for (i = 0; i < sio_add.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	std::cout << sio_add[i] << " ";
	} std::cout << " \n";
	for (i = 0; i < sap_add.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	std::cout << sap_add[i] << " " << sap_add[i + 1] << " " << sap_add[i + 2] << " ; ";
	} std::cout << " \n"; */

	// jiny pristup feasibility:
	cells = sr_add;
	cells.push_back(ijk_add);
	cells.push_back(q_add);

	//	1) Feasibility of the neighborhood pr+sr:
	if (feasibility(con, cells, alfa, beta, B)) { /*std::cout << "OK. \n";*/ }
	else { erase(ijk_add, q_add, &con); con.put(del, x_del, y_del, z_del, rad_del); return 0; }

	//	2)Feasibility of the whole configuration:
	//	if (feasibility(con, alfa, beta, B)) { /*std::cout << "OK. \n";*/ }
	//	else { erase(ijk_add, q_add, &con); con.put(del, x_del, y_del, z_del, rad_del); return 0; }

	//	3) Feasibility on the subwindow:
	//if (feasibility_subwindows(con, alfa, beta, B, nx, ny, nz, rad_del)) { /*std::cout << "OK. \n";*/ }
	//else { erase(ijk_add, q_add, &con); con.put(del, x_del, y_del, z_del, rad_del); return 0; }

	part1 = 0;
	part2 = 0;
	part3 = 0;
	l = 0;

	// 4.12.2017: domnenka, ze je nutne proverit feasibilitu i pro celou strukturu sr_add, a ne jenom pro sr_del - provede se spolecne s vypoctem part2
	// 24.10.2017: domnenka, ze i kdyz jsou neigh_add a neigh_del identicke, prece jen muze ve 3D dojit ke zmene sousedske struktury mezi nimi (???) 
	//				tedy feasibilita pro sr_del se musi overit vzdy a (ne)lze tak provest zaroven s vypoctem part2 a part3
	//if (k < neigh_add.size()) {   // nejsou-li neigh_add, neigh_del identicke musime overit feasibilitu i pro smazanou castici

	for (i = 0; i < k; i++) {				// {3} loop over the secondary particles
		ijk_i = sr_add[2 * i]; q_i = sr_add[2 * i + 1];
		con.compute_cell(d, ijk_i, q_i);  // {3} compute the cell of this primary neighbor

		if (sio_add[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap_add[3 * (sio_add[i] - 1)];
			y = sap_add[3 * (sio_add[i] - 1) + 1];
			z = sap_add[3 * (sio_add[i] - 1) + 2];
		}

		part1 = part1 + V2(c_add, d, nx, ny, nz, x, y, z);
	}

	k2 = 0;

	for (i = 0; i < sr_add.size() / 2; i++) {								// {2} loop over secondary particles 
		ijk_i = sr_add[2 * i]; q_i = sr_add[2 * i + 1];
		con.compute_cell(c, ijk_i, q_i);				// {2} compute cell of i-th secondary particle


		if (sio_add[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap_add[3 * (sio_add[i] - 1)];
			y = sap_add[3 * (sio_add[i] - 1) + 1];
			z = sap_add[3 * (sio_add[i] - 1) + 2];
		}

		for (j = i + 1; j < sr_add.size() / 2; j++) {                   // {2} loop over secondary particles with "higher order" (to avoid double counting) 
			if (are_neighbors(c, sr_add[2 * j], sr_add[2 * j + 1], &con)) {	// {2} sousedi-li bunky i a j ... (i.e. if c_i ~ c_j)
				ijk_j = sr_add[2 * j]; q_j = sr_add[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		// {2} compute cell of the j-th secondary particle

				if (sio_add[j] == 0) {
					xn = con.p[ijk_j][4 * q_j];					// computes true coordinates of j-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = sap_add[3 * (sio_add[j] - 1)];
					yn = sap_add[3 * (sio_add[j] - 1) + 1];
					zn = sap_add[3 * (sio_add[j] - 1) + 2];
				}

				part2 = part2 + V2(c, d, x, y, z, xn, yn, zn);

			} // END if (are_neighbors)
		} // END for (j; loop over secondary particles with higher order)
		c.neighbors(neigh);									           // {3} compute its neighbors ("secondary")
																	   //c.face_vertices(vert);					// {3PER} computes list of vertices for each face of the cell
																	   //fng = 0;								// {3PER}


		for (j = 0; j < neigh.size(); j++) {				           // {3} loop over neighbors of i-th secondary particle
			if (terciary(neigh[j], del, sr_add, &con)) {				 // verify if the neighbor is terciary or not (tj. nepatri mezi sekundarni castice a zaroven se take nejedna o nove pridanou castici
				find_pos(ijk, q, neigh[j], &con);			             // find its position

				tr.push_back(ijk); tr.push_back(q);			             // store this secondary neighbor to the vector of terciary particles 

				xn = con.p[ijk][4 * q]; yn = con.p[ijk][4 * q + 1]; zn = con.p[ijk][4 * q + 2];
				real_coo(x, y, z, xn, yn, zn);
				// dale jiz ale mohu pouzit skutecne souradnice bodu bunky (x,y,z):
				xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;// {3PER} urcim skutecne souradnice tohoto souseda za stenou

																  //					if (xn == con.p[ijk][3 * q] && yn == con.p[ijk][3 * q + 1] && zn == con.p[ijk][3 * q + 2]) { // {3PER} nedoslo k preklopeni terciarni castice (terciarni castice je "in")
				if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {3PER} nedoslo k preklopeni
					tio.push_back(0);
					con.compute_cell(d, ijk, q);
					part3 = part3 + V2(c, d, x, y, z, xn, yn, zn);

				}
				else {								// terciarni castice je preklopena
					k2++; tio.push_back(k2);
					tap.push_back(xn); tap.push_back(yn); tap.push_back(zn); // storing real position of this neighbor
					con.compute_cell(d, ijk, q);
					part3 = part3 + V2(c, d, x, y, z, xn, yn, zn);

				} // END if..else (terciary particle is IN vs OUT)
				l++;
			} // END if (terciary)
			  //fng = fng + vert[fng] + 1;
		} // END for (j; loop over neigbors of i-th secondary particle)
		ntr.push_back(l);
	} // END for (i; loop over secondary particles)


	  //	std::cout << "part2: " << part2 << " \n";  /////////////////////////////////////////////////////////////////////
	  //	std::cout << "part3: " << part3 << " \n";  /////////////////////////////////////////////////////////////////////

	  /* for (i = 0; i < ntr.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << ntr[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < tr.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << tr[2 * i] << " " << tr[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < tio.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << tio[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < tap.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << tap[i] << " " << tap[i + 1] << " " << tap[i + 2] << " ; ";
	  } std::cout << " \n"; */

	erase(ijk_add, q_add, &con);			// delete new particle
	con.put(del, x_del, y_del, z_del, rad_del);		// add back old particle

	find_pos(ijk_del, q_del, del, &con);	// actualize position of old particle (ijk by se zmenit nemelo, ale q ano - optimalizovat!
	con.compute_cell(c_del, ijk_del, q_del);// compute the cell of deleted particle

	part4 = 0;

	for (i = 0; i < sr_del.size() / 2; i++) {					// {4} loop over secondary particles of DEL
		ijk_i = sr_del[2 * i]; q_i = sr_del[2 * i + 1];
		con.compute_cell(c, ijk_i, q_i);  // {4} compute the cell of this primary neighbor

		if (sio_del[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap_del[3 * (sio_del[i] - 1)];
			y = sap_del[3 * (sio_del[i] - 1) + 1];
			z = sap_del[3 * (sio_del[i] - 1) + 2];
		}

		part4 = part4 + V2(c_add, d, x_del, y_del, z_del, x, y, z);

	}
	//	std::cout << "part4: " << part4 << " \n";  ////////////////////////////////////////////////////////////////////////

	// std::cout << "k1_add, k1_del, k2: " << k1_add << " " << k1_del << " " << k2 << " \n"; //////////////////////////////////

	part5 = 0;
	part6 = 0;
	for (i = 0; i < sr_add.size() / 2; i++) {                              // {5,6} loop over secondary particles 
		ijk_i = sr_add[2 * i]; q_i = sr_add[2 * i + 1];
		if (sio_add[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap_add[3 * (sio_add[i] - 1)];
			y = sap_add[3 * (sio_add[i] - 1) + 1];
			z = sap_add[3 * (sio_add[i] - 1) + 2];
		}
		con.compute_cell(c, ijk_i, q_i);			// {5,6} compute the cell of i-th primary neighbor

		for (j = i + 1; j < sr_add.size() / 2; j++) {				  // {5} loop over secondary particles with "higher order" (to avoid double counting)  	
			if (are_neighbors(c, sr_add[2 * j], sr_add[2 * j + 1], &con)) {	  // if i and j are neighbors... (i.e. if c_i ~ c_j)
				ijk_j = sr_add[2 * j]; q_j = sr_add[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		  // compute cell j

				if (sio_add[j] == 0) {
					xn = con.p[ijk_j][4 * q_j];					// computes coordinates of j-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = sap_add[3 * (sio_add[j] - 1)];
					yn = sap_add[3 * (sio_add[j] - 1) + 1];
					zn = sap_add[3 * (sio_add[j] - 1) + 2];
				}

				part5 = part5 + V2(c, d, x, y, z, xn, yn, zn);
			}
		}

		if (i == 0) {
			for (j = 0; j < ntr[0]; j++) {							// {6} i=0; loop over its terciary particles		
				ijk_j = tr[2 * j]; q_j = tr[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		  // compute cell of this terciary particle

				if (tio[j] == 0) {
					xn = con.p[ijk_j][4 * q_j];								// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = tap[3 * (tio[j] - 1)];
					yn = tap[3 * (tio[j] - 1) + 1];
					zn = tap[3 * (tio[j] - 1) + 2];
				}

				part6 = part6 + V2(c, d, x, y, z, xn, yn, zn);
			}
		}
		else {
			for (j = 0; j < ntr[i] - ntr[i - 1]; j++) {				// {6} i>0; loop over its terciary particles
				k = ntr[i - 1] + j;
				ijk_j = tr[2 * k]; q_j = tr[2 * k + 1];
				con.compute_cell(d, ijk_j, q_j);

				if (tio[k] == 0) {
					xn = con.p[ijk_j][4 * q_j];								// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = tap[3 * (tio[k] - 1)];
					yn = tap[3 * (tio[k] - 1) + 1];
					zn = tap[3 * (tio[k] - 1) + 2];
				}

				part6 = part6 + V2(c, d, x, y, z, xn, yn, zn);

			}
		} // END if..else (i=0 vs i>0)

	} // END for (i; loop over secondary particles)

	  //	std::cout << "part5: " << part5 << " \n";  /////////////////////////////////////////////////////////////////////
	  //	std::cout << "part6: " << part6 << " \n";  /////////////////////////////////////////////////////////////////////

	  //	std::cout << "M: part1: " << -part1 << " | part2: " << -part2 << " | part3: " << -part3 << " | part4: " << part4 << " | part5: " << part5 << " | part6: " << part6 << " \n";

	return exp(theta*(-part1 - part2 - part3 + part4 + part5 + part6));
}





// 770 - 1227

double try_move(int &ijk_del, int &q_del, double &nx, double &ny, double &nz, container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, const double &iota)
{
	// [in]		ijk_del,q_del	position of the particle to be deleted.
	// [in]		x,y,z			the ID number and coordinates of the particle to be added.  -- misto tohoto by stacilo predat hodnotu sigma
	// [in]		con				the container with stored particles.
	// [in]		theta			parameter
	// [in]		alfa,beta,B		hardcore parameters.

	unsigned int i, j, k, l, k1_add, k1_del, k2, fng;
	int ijk, q, ijk_add, q_add, del, ijk_i, q_i, ijk_j, q_j;
	double part1, part2, part3, part4, part5, part6;
	double x_del, y_del, z_del, x, y, z, xn, yn, zn;
	double vol, rad_del;

	std::vector<int> sr_add, sr_del, tr; // misto particle_order class snadnejsi pouzivat vektory pro ukladani ijk,q castic
	std::vector<int> sio_add, sio_del, tio;
	std::vector<double> sap_add, sap_del, tap;
	std::vector<int> ntr, cells;
	std::vector<int> neigh_del, neigh_add, neigh, neighi;  // vektor pro prirazeni ID sousedu
	std::vector<int> vert, verti;
	voronoicell_neighbor c_del, c_add, c, d;  // bunka s informaci o sousedech


	x_del = con.p[ijk_del][4 * q_del];      // coordinates of deleted particle
	y_del = con.p[ijk_del][4 * q_del + 1];
	z_del = con.p[ijk_del][4 * q_del + 2];
	rad_del = con.p[ijk_del][4 * q_del + 3];

	del = con.id[ijk_del][q_del];			// ID of deleted particle

	con.compute_cell(c_del, ijk_del, q_del);	// compute the cell of deleted particle
	c_del.neighbors(neigh_del);					// compute its neighbors (in following "primary")
												//n_del = neigh_del.size();					// determine number of neighbors of deleted particle
												//c_del.face_vertices(vert);					// {PER} computes list of vertices for each face of the cell
												//fng = 0;

	erase(ijk_del, q_del, &con);				// delete old particle
												// std::cout << "del \n";     //////////////////////////////////////////////////////////////////////////////////////////////

	k1_del = 0;
	for (i = 0; i < neigh_del.size(); i++) {
		find_pos(ijk, q, neigh_del[i], &con);

		sr_del.push_back(ijk); sr_del.push_back(q);
		// std::cout << con.p[ijk][4 * q] << " " << con.p[ijk][4 * q + 1] << " " << con.p[ijk][4 * q + 2] << " ; "; /////////////////

		xn = con.p[ijk][4 * q]; yn = con.p[ijk][4 * q + 1]; zn = con.p[ijk][4 * q + 2];
		real_coo(x_del, y_del, z_del, xn, yn, zn);

		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {PER} nedoslo k preklopeni
			sio_del.push_back(0);					// neighbor is "in"
		}
		else {
			k1_del++;
			sio_del.push_back(k1_del);					// neighbor is k-th "out"
			sap_del.push_back(xn); sap_del.push_back(yn); sap_del.push_back(zn);  // storing real position of this neighbor
		}

		//fng = fng + vert[fng] + 1;
	} // END (i; loop over neighbors of deleted particle)
	  //vert.clear();

	  // std::cout << "\n";     //////////////////////////////////////////////////////////////////////////////////////////////
	  /* for (i = 0; i < neigh_del.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << neigh_del[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sr_del.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << sr_del[2 * i] << " " << sr_del[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < sio_del.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sio_del[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sap_del.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sap_del[i] << " " << sap_del[i + 1] << " " << sap_del[i + 2] << " ; ";
	  } std::cout << " \n"; */

	con.put(del, nx, ny, nz, rad_del);					// add new particle
														// co kdyz je nektera z nx, ny, nz mimo okno? put si asi diky periodicite poradi, ale poradi si i zbytek???
														//   konkretne poradi si fce face_dist? ne fci face_dist se musi predat pozice v oknë!!!!!!!!!!!
														// je potreba premapovat nx,ny,nz do okna ---> implementovano jiz ve fci bdma_step

	find_pos(ijk_add, q_add, del, &con);			// determine position of new particle

	if (overlap_f(con, iota)) {}
	else { erase(ijk_add, q_add, &con); con.put(del, x_del, y_del, z_del, rad_del); return 0; }


	con.compute_cell(c_add, ijk_add, q_add);	// compute cell for new added particle
	c_add.neighbors(neigh_add);					// compute its neighbors	
												//n_add = neigh_add.size();					// determine number of neighbors of added particle      !!!mozna zbytecne
												//c_add.face_vertices(vert);					// {F,1PER} computes list of vertices for each face of the cell
												//fng = 0;								// {F,1PER}

												// std::cout << ijk_add << " " << q_add << "\n"; ///////////////////////////////////////////////////////////////////////////

	k1_add = 0;

	vol = c_add.volume();
	part1 = 0;

	for (i = 0; i < neigh_add.size(); i++) {	// loop over the neighbors
		find_pos(ijk, q, neigh_add[i], &con);	// find position of the i-th neighbor

												// std::cout << con.p[ijk][4 * q] << " " << con.p[ijk][4 * q + 1] << " " << con.p[ijk][4 * q + 2] << " ; "; /////////////////

		sr_add.push_back(ijk); sr_add.push_back(q); // add this neighbor to the second particle vector

		xn = con.p[ijk][4 * q]; yn = con.p[ijk][4 * q + 1]; zn = con.p[ijk][4 * q + 2];
		real_coo(nx, ny, nz, xn, yn, zn);

		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {PER} nedoslo k preklopeni
																			//part1 = part1 + V2(c_add, d);			// {1} add the value of V2 function to the first sum
			sio_add.push_back(0);					// neighbor is "in"
		}
		else {
			//part1 = part1 + V2(c_add, d, nx, ny, nz, xn, yn, zn);
			k1_add++;
			sio_add.push_back(k1_add);					// neighbor is k-th "out"
			sap_add.push_back(xn); sap_add.push_back(yn); sap_add.push_back(zn);  // storing real position of this neighbor
		}
		//fng = fng + vert[fng] + 1;			// {1PER} set actual position in vector of face vertices
	} // END for (i; loop over neighbors of added particle)

	  // std::cout << "\n";    ////////////////////////////////////////////////////////////////////////////////////////////////////

	  /* std::cout << "add \n";     //////////////////////////////////////////////////////////////////////////////////////////////
	  for (i = 0; i < neigh_add.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << neigh_add[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sr_add.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << sr_add[2 * i] << " " << sr_add[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < sio_add.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sio_add[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sap_add.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sap_add[i] << " " << sap_add[i + 1] << " " << sap_add[i + 2] << " ; ";
	  } std::cout << " \n"; */

	  //	std::cout << "part1: " << part1 << " \n";  ////////////////////////////////////////////////////////////////////////////////

	  // B - merge datovych struktur
	k = neigh_add.size();
	merge(neigh_add, neigh_del, sr_add, sr_del, k1_add, sio_add, sio_del, sap_add, sap_del);
	// slouceni se provede do vektoru sr_add, sio_add, sap_add
	// vektory sr_del, sio_del, sap_del se mi jeste budou hodit (jeste nebyla overena pripustnost pro vsechny zmenene castice)
	// vektor neigh_del byl podmnozinou vektoru neigh_add, pokud stale plati k == neigh_add.size() i po slouceni (tj. nebylo nic pridano)

	/* std::cout << "merge \n"; /////////////////////////////////////////////////////////////////////////////////////////////////
	for (i = 0; i < sr_add.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	std::cout << sr_add[2 * i] << " " << sr_add[2 * i + 1] << " ; ";
	} std::cout << " \n";
	for (i = 0; i < sio_add.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	std::cout << sio_add[i] << " ";
	} std::cout << " \n";
	for (i = 0; i < sap_add.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	std::cout << sap_add[i] << " " << sap_add[i + 1] << " " << sap_add[i + 2] << " ; ";
	} std::cout << " \n"; */

	// jiny pristup feasibility:
	cells = sr_add;
	cells.push_back(ijk_add);
	cells.push_back(q_add);

	//	1) Feasibility of the neighborhood pr+sr:
		if (feasibility(con, cells, alfa, beta, B)) { /*std::cout << "OK. \n";*/ }
		else { erase(ijk_add, q_add, &con); con.put(del, x_del, y_del, z_del, rad_del); return 0; }

	//	2)Feasibility of the whole configuration:
	//	if (feasibility(con, alfa, beta, B)) { /*std::cout << "OK. \n";*/ }
	//	else { erase(ijk_add, q_add, &con); con.put(del, x_del, y_del, z_del, rad_del); return 0; }

	//	3) Feasibility on the subwindow:
	//if (feasibility_subwindows(con, alfa, beta, B, nx, ny, nz, rad_del)) { /*std::cout << "OK. \n";*/ }
	//else { erase(ijk_add, q_add, &con); con.put(del, x_del, y_del, z_del, rad_del); return 0; }

	part1 = 0;
	part2 = 0;
	part3 = 0;
	l = 0;

	// 4.12.2017: domnenka, ze je nutne proverit feasibilitu i pro celou strukturu sr_add, a ne jenom pro sr_del - provede se spolecne s vypoctem part2
	// 24.10.2017: domnenka, ze i kdyz jsou neigh_add a neigh_del identicke, prece jen muze ve 3D dojit ke zmene sousedske struktury mezi nimi (???) 
	//				tedy feasibilita pro sr_del se musi overit vzdy a (ne)lze tak provest zaroven s vypoctem part2 a part3
	//if (k < neigh_add.size()) {   // nejsou-li neigh_add, neigh_del identicke musime overit feasibilitu i pro smazanou castici

	for (i = 0; i < k; i++) {				// {3} loop over the secondary particles
		ijk_i = sr_add[2 * i]; q_i = sr_add[2 * i + 1];
		con.compute_cell(d, ijk_i, q_i);  // {3} compute the cell of this primary neighbor

		if (sio_add[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap_add[3 * (sio_add[i] - 1)];
			y = sap_add[3 * (sio_add[i] - 1) + 1];
			z = sap_add[3 * (sio_add[i] - 1) + 2];
		}

		part1 = part1 + V2(c_add, d, nx, ny, nz, x, y, z);
	}

	k2 = 0;

	for (i = 0; i < sr_add.size() / 2; i++) {								// {2} loop over secondary particles 
		ijk_i = sr_add[2 * i]; q_i = sr_add[2 * i + 1];
		con.compute_cell(c, ijk_i, q_i);				// {2} compute cell of i-th secondary particle


		if (sio_add[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap_add[3 * (sio_add[i] - 1)];
			y = sap_add[3 * (sio_add[i] - 1) + 1];
			z = sap_add[3 * (sio_add[i] - 1) + 2];
		}

		for (j = i + 1; j < sr_add.size() / 2; j++) {                   // {2} loop over secondary particles with "higher order" (to avoid double counting) 
			if (are_neighbors(c, sr_add[2 * j], sr_add[2 * j + 1], &con)) {	// {2} sousedi-li bunky i a j ... (i.e. if c_i ~ c_j)
				ijk_j = sr_add[2 * j]; q_j = sr_add[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		// {2} compute cell of the j-th secondary particle

				if (sio_add[j] == 0) {
					xn = con.p[ijk_j][4 * q_j];					// computes true coordinates of j-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = sap_add[3 * (sio_add[j] - 1)];
					yn = sap_add[3 * (sio_add[j] - 1) + 1];
					zn = sap_add[3 * (sio_add[j] - 1) + 2];
				}

				part2 = part2 + V2(c, d, x, y, z, xn, yn, zn);

			} // END if (are_neighbors)
		} // END for (j; loop over secondary particles with higher order)
		c.neighbors(neigh);									           // {3} compute its neighbors ("secondary")
																	   //c.face_vertices(vert);					// {3PER} computes list of vertices for each face of the cell
																	   //fng = 0;								// {3PER}


		for (j = 0; j < neigh.size(); j++) {				           // {3} loop over neighbors of i-th secondary particle
			if (terciary(neigh[j], del, sr_add, &con)) {				 // verify if the neighbor is terciary or not (tj. nepatri mezi sekundarni castice a zaroven se take nejedna o nove pridanou castici
				find_pos(ijk, q, neigh[j], &con);			             // find its position

				tr.push_back(ijk); tr.push_back(q);			             // store this secondary neighbor to the vector of terciary particles 

				xn = con.p[ijk][4 * q]; yn = con.p[ijk][4 * q + 1]; zn = con.p[ijk][4 * q + 2];
				real_coo(x, y, z, xn, yn, zn);
				// dale jiz ale mohu pouzit skutecne souradnice bodu bunky (x,y,z):
				xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;// {3PER} urcim skutecne souradnice tohoto souseda za stenou

																  //					if (xn == con.p[ijk][3 * q] && yn == con.p[ijk][3 * q + 1] && zn == con.p[ijk][3 * q + 2]) { // {3PER} nedoslo k preklopeni terciarni castice (terciarni castice je "in")
				if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {3PER} nedoslo k preklopeni
					tio.push_back(0);
					con.compute_cell(d, ijk, q);
					part3 = part3 + V2(c, d, x, y, z, xn, yn, zn);

				}
				else {								// terciarni castice je preklopena
					k2++; tio.push_back(k2);
					tap.push_back(xn); tap.push_back(yn); tap.push_back(zn); // storing real position of this neighbor
					con.compute_cell(d, ijk, q);
					part3 = part3 + V2(c, d, x, y, z, xn, yn, zn);

				} // END if..else (terciary particle is IN vs OUT)
				l++;
			} // END if (terciary)
			  //fng = fng + vert[fng] + 1;
		} // END for (j; loop over neigbors of i-th secondary particle)
		ntr.push_back(l);
	} // END for (i; loop over secondary particles)


	  //	std::cout << "part2: " << part2 << " \n";  /////////////////////////////////////////////////////////////////////
	  //	std::cout << "part3: " << part3 << " \n";  /////////////////////////////////////////////////////////////////////

	  /* for (i = 0; i < ntr.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << ntr[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < tr.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << tr[2 * i] << " " << tr[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < tio.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << tio[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < tap.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << tap[i] << " " << tap[i + 1] << " " << tap[i + 2] << " ; ";
	  } std::cout << " \n"; */

	erase(ijk_add, q_add, &con);			// delete new particle
	con.put(del, x_del, y_del, z_del, rad_del);		// add back old particle

	find_pos(ijk_del, q_del, del, &con);	// actualize position of old particle (ijk by se zmenit nemelo, ale q ano - optimalizovat!
	con.compute_cell(c_del, ijk_del, q_del);// compute the cell of deleted particle

	part4 = 0;

	for (i = 0; i < sr_del.size() / 2; i++) {					// {4} loop over secondary particles of DEL
		ijk_i = sr_del[2 * i]; q_i = sr_del[2 * i + 1];
		con.compute_cell(c, ijk_i, q_i);  // {4} compute the cell of this primary neighbor

		if (sio_del[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap_del[3 * (sio_del[i] - 1)];
			y = sap_del[3 * (sio_del[i] - 1) + 1];
			z = sap_del[3 * (sio_del[i] - 1) + 2];
		}

		part4 = part4 + V2(c_add, d, x_del, y_del, z_del, x, y, z);

	}
	//	std::cout << "part4: " << part4 << " \n";  ////////////////////////////////////////////////////////////////////////

	// std::cout << "k1_add, k1_del, k2: " << k1_add << " " << k1_del << " " << k2 << " \n"; //////////////////////////////////

	part5 = 0;
	part6 = 0;
	for (i = 0; i < sr_add.size() / 2; i++) {                              // {5,6} loop over secondary particles 
		ijk_i = sr_add[2 * i]; q_i = sr_add[2 * i + 1];
		if (sio_add[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap_add[3 * (sio_add[i] - 1)];
			y = sap_add[3 * (sio_add[i] - 1) + 1];
			z = sap_add[3 * (sio_add[i] - 1) + 2];
		}
		con.compute_cell(c, ijk_i, q_i);			// {5,6} compute the cell of i-th primary neighbor

		for (j = i + 1; j < sr_add.size() / 2; j++) {				  // {5} loop over secondary particles with "higher order" (to avoid double counting)  	
			if (are_neighbors(c, sr_add[2 * j], sr_add[2 * j + 1], &con)) {	  // if i and j are neighbors... (i.e. if c_i ~ c_j)
				ijk_j = sr_add[2 * j]; q_j = sr_add[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		  // compute cell j

				if (sio_add[j] == 0) {
					xn = con.p[ijk_j][4 * q_j];					// computes coordinates of j-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = sap_add[3 * (sio_add[j] - 1)];
					yn = sap_add[3 * (sio_add[j] - 1) + 1];
					zn = sap_add[3 * (sio_add[j] - 1) + 2];
				}

				part5 = part5 + V2(c, d, x, y, z, xn, yn, zn);
			}
		}

		if (i == 0) {
			for (j = 0; j < ntr[0]; j++) {							// {6} i=0; loop over its terciary particles		
				ijk_j = tr[2 * j]; q_j = tr[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		  // compute cell of this terciary particle

				if (tio[j] == 0) {
					xn = con.p[ijk_j][4 * q_j];								// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = tap[3 * (tio[j] - 1)];
					yn = tap[3 * (tio[j] - 1) + 1];
					zn = tap[3 * (tio[j] - 1) + 2];
				}

				part6 = part6 + V2(c, d, x, y, z, xn, yn, zn);
			}
		}
		else {
			for (j = 0; j < ntr[i] - ntr[i - 1]; j++) {				// {6} i>0; loop over its terciary particles
				k = ntr[i - 1] + j;
				ijk_j = tr[2 * k]; q_j = tr[2 * k + 1];
				con.compute_cell(d, ijk_j, q_j);

				if (tio[k] == 0) {
					xn = con.p[ijk_j][4 * q_j];								// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = tap[3 * (tio[k] - 1)];
					yn = tap[3 * (tio[k] - 1) + 1];
					zn = tap[3 * (tio[k] - 1) + 2];
				}

				part6 = part6 + V2(c, d, x, y, z, xn, yn, zn);

			}
		} // END if..else (i=0 vs i>0)

	} // END for (i; loop over secondary particles)

	  //	std::cout << "part5: " << part5 << " \n";  /////////////////////////////////////////////////////////////////////
	  //	std::cout << "part6: " << part6 << " \n";  /////////////////////////////////////////////////////////////////////

	  //	std::cout << "M: part1: " << -part1 << " | part2: " << -part2 << " | part3: " << -part3 << " | part4: " << part4 << " | part5: " << part5 << " | part6: " << part6 << " \n";

	return exp(theta*(-part1 - part2 - part3 + part4 + part5 + part6));
}

// zmena polomeru odpovida v podstate kroku move, prtz vysledkem zmeny polomeru muze byt i zmena sousedicich castic
double try_change_rad(int &ijk_chr, int &q_chr, double &nrad, container_poly &con, const double &theta, const double &alfa, const double &beta, const double &B, const double &iota)
{
	// [in]		ijk_chr,q_chr	position of the particle.
	// [in]		rad				suggestion of the new particle radius
	// [in]		con				the container with stored particles.
	// [in]		theta			parameter
	// [in]		alfa,beta,B		hardcore parameters.

	double radc;
	unsigned int i, j, k, l, k1_add, k1_del, k2, fng;
	int ijk, q, chr, ijk_i, q_i, ijk_j, q_j;
	double part1, part2, part3, part4, part5, part6;
	double xc, yc, zc, x, y, z, xn, yn, zn;

	std::vector<int> sr_b, sr_a, tr; // misto particle_order class snadnejsi pouzivat vektory pro ukladani ijk,q castic
	std::vector<int> sio_b, sio_a, tio;
	std::vector<double> sap_b, sap_a, tap;
	std::vector<int> ntr, cells;
	std::vector<int> neigh_b, neigh_a, neighi;  // vektor pro prirazeni ID sousedu
	std::vector<int> vert, verti;
	voronoicell_neighbor c_b, c_a, c, d;  // bunka s informaci o sousedech

	xc = con.p[ijk_chr][4 * q_chr];      // coordinates of changed particle
	yc = con.p[ijk_chr][4 * q_chr + 1];
	zc = con.p[ijk_chr][4 * q_chr + 2];
	radc = con.p[ijk_chr][4 * q_chr + 3];	// radius of changed particle
	chr = con.id[ijk_chr][q_chr];			// ID of changed particle

											// before change:
	con.compute_cell(c_b, ijk_chr, q_chr);	// compute the cell of changed particle
	c_b.neighbors(neigh_b);						// compute its neighbors (in following "primary")
												//n_del = neigh_del.size();					// determine number of neighbors of particle
												//c_b.face_vertices(vert);					// {PER} computes list of vertices for each face of the cell
												//fng = 0;

	k1_del = 0;
	for (i = 0; i < neigh_b.size(); i++) {
		find_pos(ijk, q, neigh_b[i], &con);

		sr_b.push_back(ijk); sr_b.push_back(q);
		// std::cout << con.p[ijk][4 * q] << " " << con.p[ijk][4 * q + 1] << " " << con.p[ijk][4 * q + 2] << " ; "; /////////////////

		xn = con.p[ijk][4 * q]; yn = con.p[ijk][4 * q + 1]; zn = con.p[ijk][4 * q + 2];
		real_coo(xc, yc, zc, xn, yn, zn);

		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {PER} nedoslo k preklopeni
			sio_b.push_back(0);					// neighbor is "in"
		}
		else {
			k1_del++;
			sio_b.push_back(k1_del);					// neighbor is k-th "out"
			sap_b.push_back(xn); sap_b.push_back(yn); sap_b.push_back(zn);  // storing real position of this neighbor
		}

		//fng = fng + vert[fng] + 1;
	} // END (i; loop over neighbors of deleted particle)
	  //vert.clear();

	  // std::cout << "\n";     //////////////////////////////////////////////////////////////////////////////////////////////
	  /* for (i = 0; i < neigh_del.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << neigh_del[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sr_del.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << sr_del[2 * i] << " " << sr_del[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < sio_del.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sio_del[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sap_del.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sap_del[i] << " " << sap_del[i + 1] << " " << sap_del[i + 2] << " ; ";
	  } std::cout << " \n"; */

	con.p[ijk_chr][4 * q_chr + 3] = nrad;	// zmena polomeru ------------------------------------------------------------------------------------------------


	if (overlap_f(con,iota)) {}
	else { con.p[ijk_chr][4 * q_chr + 3] = radc; return 0; }

											// after change:
	con.compute_cell(c_a, ijk_chr, q_chr);	// compute cell for new added particle
	c_a.neighbors(neigh_a);					// compute its neighbors	
											//n_add = neigh_add.size();					// determine number of neighbors of added particle      !!!mozna zbytecne
											//c_a.face_vertices(vert);					// {F,1PER} computes list of vertices for each face of the cell
											//fng = 0;								// {F,1PER}

											// std::cout << ijk_add << " " << q_add << "\n"; ///////////////////////////////////////////////////////////////////////////

	k1_add = 0;
	for (i = 0; i < neigh_a.size(); i++) {	// loop over the neighbors
		find_pos(ijk, q, neigh_a[i], &con);	// find position of the i-th neighbor

											// std::cout << con.p[ijk][4 * q] << " " << con.p[ijk][4 * q + 1] << " " << con.p[ijk][4 * q + 2] << " ; "; /////////////////

		sr_a.push_back(ijk); sr_a.push_back(q); // add this neighbor to the second particle vector

		xn = con.p[ijk][4 * q]; yn = con.p[ijk][4 * q + 1]; zn = con.p[ijk][4 * q + 2];
		real_coo(xc, yc, zc, xn, yn, zn);

		if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {PER} nedoslo k preklopeni
																			//part1 = part1 + V2(c_add, d);			// {1} add the value of V2 function to the first sum
			sio_a.push_back(0);					// neighbor is "in"
		}
		else {
			//part1 = part1 + V2(c_add, d, nx, ny, nz, xn, yn, zn);
			k1_add++;
			sio_a.push_back(k1_add);					// neighbor is k-th "out"
			sap_a.push_back(xn); sap_a.push_back(yn); sap_a.push_back(zn);  // storing real position of this neighbor
		}
		//fng = fng + vert[fng] + 1;			// {1PER} set actual position in vector of face vertices
	} // END for (i; loop over neighbors of added particle)

	  // std::cout << "\n";    ////////////////////////////////////////////////////////////////////////////////////////////////////

	  /* std::cout << "add \n";     //////////////////////////////////////////////////////////////////////////////////////////////
	  for (i = 0; i < neigh_add.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << neigh_add[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sr_add.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << sr_add[2 * i] << " " << sr_add[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < sio_add.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sio_add[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < sap_add.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << sap_add[i] << " " << sap_add[i + 1] << " " << sap_add[i + 2] << " ; ";
	  } std::cout << " \n"; */

	  //	std::cout << "part1: " << part1 << " \n";  ////////////////////////////////////////////////////////////////////////////////

	  // B - merge datovych struktur
	k = neigh_a.size();
	merge(neigh_a, neigh_b, sr_a, sr_b, k1_add, sio_a, sio_b, sap_a, sap_b);
	// slouceni se provede do vektoru sr_add, sio_add, sap_add
	// vektory sr_del, sio_del, sap_del se mi jeste budou hodit (jeste nebyla overena pripustnost pro vsechny zmenene castice)
	// vektor neigh_del byl podmnozinou vektoru neigh_add, pokud stale plati k == neigh_add.size() i po slouceni (tj. nebylo nic pridano)

	/* std::cout << "merge \n"; /////////////////////////////////////////////////////////////////////////////////////////////////
	for (i = 0; i < sr_add.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	std::cout << sr_add[2 * i] << " " << sr_add[2 * i + 1] << " ; ";
	} std::cout << " \n";
	for (i = 0; i < sio_add.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	std::cout << sio_add[i] << " ";
	} std::cout << " \n";
	for (i = 0; i < sap_add.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	std::cout << sap_add[i] << " " << sap_add[i + 1] << " " << sap_add[i + 2] << " ; ";
	} std::cout << " \n"; */

	// jiny pristup feasibility:
	cells = sr_a;
	cells.push_back(ijk_chr);
	cells.push_back(q_chr);

	//	1) Feasibility of the neighborhood pr+sr:
		if (feasibility(con, cells, alfa, beta, B)) { /*std::cout << "OK. \n";*/ }
		else { con.p[ijk_chr][4 * q_chr + 3] = radc; return 0; }

	//	2) Feasibility of the whole configuration:
	//	if (feasibility(con, alfa, beta, B)) { /*std::cout << "OK. \n";*/ }
	//	else { con.p[ijk_chr][4 * q_chr + 3] = radc; return 0; }

	//	3) Feasibility on the subwindow:
	//if (feasibility_subwindows(con, alfa, beta, B, xc, yc, zc, nrad)) { /*std::cout << "OK. \n";*/ }
	//else { con.p[ijk_chr][4 * q_chr + 3] = radc; return 0; }


	part1 = 0;
	part2 = 0;
	part3 = 0;

	for (i = 0; i < k; i++) {				// {3} loop over the secondary particles
		ijk_i = sr_a[2 * i]; q_i = sr_a[2 * i + 1];
		con.compute_cell(d, ijk_i, q_i);  // {3} compute the cell of this primary neighbor

		if (sio_a[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap_a[3 * (sio_a[i] - 1)];
			y = sap_a[3 * (sio_a[i] - 1) + 1];
			z = sap_a[3 * (sio_a[i] - 1) + 2];
		}

		part1 = part1 + V2(c_a, d, xc, yc, zc, x, y, z);
	}

	k2 = 0;
	l = 0;

	for (i = 0; i < sr_a.size() / 2; i++) {								// {2} loop over secondary particles 
		ijk_i = sr_a[2 * i]; q_i = sr_a[2 * i + 1];
		con.compute_cell(c, ijk_i, q_i);				// {2} compute cell of i-th secondary particle


		if (sio_a[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap_a[3 * (sio_a[i] - 1)];
			y = sap_a[3 * (sio_a[i] - 1) + 1];
			z = sap_a[3 * (sio_a[i] - 1) + 2];
		}

		for (j = i + 1; j < sr_a.size() / 2; j++) {                   // {2} loop over secondary particles with "higher order" (to avoid double counting) 
			if (are_neighbors(c, sr_a[2 * j], sr_a[2 * j + 1], &con)) {	// {2} sousedi-li bunky i a j ... (i.e. if c_i ~ c_j)
				ijk_j = sr_a[2 * j]; q_j = sr_a[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		// {2} compute cell of the j-th secondary particle

				if (sio_a[j] == 0) {
					xn = con.p[ijk_j][4 * q_j];					// computes true coordinates of j-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = sap_a[3 * (sio_a[j] - 1)];
					yn = sap_a[3 * (sio_a[j] - 1) + 1];
					zn = sap_a[3 * (sio_a[j] - 1) + 2];
				}

				part2 = part2 + V2(c, d, x, y, z, xn, yn, zn);

			} // END if (are_neighbors)
		} // END for (j; loop over secondary particles with higher order)
		c.neighbors(neighi);									           // {3} compute its neighbors ("secondary")
																		   //c.face_vertices(verti);					// {3PER} computes list of vertices for each face of the cell
																		   //fng = 0;								// {3PER}

		for (j = 0; j < neighi.size(); j++) {				           // {3} loop over neighbors of i-th secondary particle
			if (terciary(neighi[j], chr, sr_a, &con)) {				 // verify if the neighbor is terciary or not (tj. nepatri mezi sekundarni castice a zaroven se take nejedna o nove pridanou castici
				find_pos(ijk, q, neighi[j], &con);			             // find its position

				tr.push_back(ijk); tr.push_back(q);			             // store this secondary neighbor to the vector of terciary particles 

				xn = con.p[ijk][4 * q]; yn = con.p[ijk][4 * q + 1]; zn = con.p[ijk][4 * q + 2];
				real_coo(x, y, z, xn, yn, zn);

				if (xn > 0 && xn < 1 && yn > 0 && yn < 1 && zn > 0 && zn < 1) {		// {3PER} nedoslo k preklopeni
					tio.push_back(0);
					con.compute_cell(d, ijk, q);
					part3 = part3 + V2(c, d, x, y, z, xn, yn, zn);

				}
				else {								// terciarni castice je preklopena
					k2++; tio.push_back(k2);
					tap.push_back(xn); tap.push_back(yn); tap.push_back(zn); // storing real position of this neighbor
					con.compute_cell(d, ijk, q);
					part3 = part3 + V2(c, d, x, y, z, xn, yn, zn);

				} // END if..else (terciary particle is IN vs OUT)
				l++;
			} // END if (terciary)
			  //fng = fng + verti[fng] + 1;
		} // END for (j; loop over neigbors of i-th secondary particle)
		ntr.push_back(l);
	} // END for (i; loop over secondary particles)


	  //	std::cout << "part2: " << part2 << " \n";  /////////////////////////////////////////////////////////////////////
	  //	std::cout << "part3: " << part3 << " \n";  /////////////////////////////////////////////////////////////////////

	  /* for (i = 0; i < ntr.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << ntr[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < tr.size() / 2; i++) {  ///////////////////////////////////////////////////////////////////////////////
	  std::cout << tr[2 * i] << " " << tr[2 * i + 1] << " ; ";
	  } std::cout << " \n";
	  for (i = 0; i < tio.size(); i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << tio[i] << " ";
	  } std::cout << " \n";
	  for (i = 0; i < tap.size() / 3; i++) {  /////////////////////////////////////////////////////////////////////////////////
	  std::cout << tap[i] << " " << tap[i + 1] << " " << tap[i + 2] << " ; ";
	  } std::cout << " \n"; */

	con.p[ijk_chr][4 * q_chr + 3] = radc;	// change radius back


	part4 = 0;

	for (i = 0; i < sr_b.size() / 2; i++) {					// {4} loop over secondary particles of DEL
		ijk_i = sr_b[2 * i]; q_i = sr_b[2 * i + 1];
		con.compute_cell(c, ijk_i, q_i);  // {4} compute the cell of this primary neighbor

		if (sio_b[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// {F} computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap_b[3 * (sio_b[i] - 1)];
			y = sap_b[3 * (sio_b[i] - 1) + 1];
			z = sap_b[3 * (sio_b[i] - 1) + 2];
		}

		part4 = part4 + V2(c_b, d, xc, yc, zc, x, y, z);

	}
	//	std::cout << "part4: " << part4 << " \n";  ////////////////////////////////////////////////////////////////////////

	// std::cout << "k1_add, k1_del, k2: " << k1_add << " " << k1_del << " " << k2 << " \n"; //////////////////////////////////

	part5 = 0;
	part6 = 0;
	for (i = 0; i < sr_a.size() / 2; i++) {                              // {5,6} loop over secondary particles 
		ijk_i = sr_a[2 * i]; q_i = sr_a[2 * i + 1];
		if (sio_a[i] == 0) {
			x = con.p[ijk_i][4 * q_i];					// computes coordinates of i-th particle
			y = con.p[ijk_i][4 * q_i + 1];
			z = con.p[ijk_i][4 * q_i + 2];
		}
		else {
			x = sap_a[3 * (sio_a[i] - 1)];
			y = sap_a[3 * (sio_a[i] - 1) + 1];
			z = sap_a[3 * (sio_a[i] - 1) + 2];
		}
		con.compute_cell(c, ijk_i, q_i);			// {5,6} compute the cell of i-th primary neighbor

		for (j = i + 1; j < sr_a.size() / 2; j++) {				  // {5} loop over secondary particles with "higher order" (to avoid double counting)  	
			if (are_neighbors(c, sr_a[2 * j], sr_a[2 * j + 1], &con)) {	  // if i and j are neighbors... (i.e. if c_i ~ c_j)
				ijk_j = sr_a[2 * j]; q_j = sr_a[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		  // compute cell j

				if (sio_a[j] == 0) {
					xn = con.p[ijk_j][4 * q_j];					// computes coordinates of j-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = sap_a[3 * (sio_a[j] - 1)];
					yn = sap_a[3 * (sio_a[j] - 1) + 1];
					zn = sap_a[3 * (sio_a[j] - 1) + 2];
				}

				part5 = part5 + V2(c, d, x, y, z, xn, yn, zn);
			}
		}

		if (i == 0) {
			for (j = 0; j < ntr[0]; j++) {							// {6} i=0; loop over its terciary particles		
				ijk_j = tr[2 * j]; q_j = tr[2 * j + 1];
				con.compute_cell(d, ijk_j, q_j);		  // compute cell of this terciary particle

				if (tio[j] == 0) {
					xn = con.p[ijk_j][4 * q_j];								// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = tap[3 * (tio[j] - 1)];
					yn = tap[3 * (tio[j] - 1) + 1];
					zn = tap[3 * (tio[j] - 1) + 2];
				}

				part6 = part6 + V2(c, d, x, y, z, xn, yn, zn);
			}
		}
		else {
			for (j = 0; j < ntr[i] - ntr[i - 1]; j++) {				// {6} i>0; loop over its terciary particles
				k = ntr[i - 1] + j;
				ijk_j = tr[2 * k]; q_j = tr[2 * k + 1];
				con.compute_cell(d, ijk_j, q_j);

				if (tio[k] == 0) {
					xn = con.p[ijk_j][4 * q_j];								// {F} computes coordinates of i-th particle
					yn = con.p[ijk_j][4 * q_j + 1];
					zn = con.p[ijk_j][4 * q_j + 2];
				}
				else {
					xn = tap[3 * (tio[k] - 1)];
					yn = tap[3 * (tio[k] - 1) + 1];
					zn = tap[3 * (tio[k] - 1) + 2];
				}

				part6 = part6 + V2(c, d, x, y, z, xn, yn, zn);

			}
		} // END if..else (i=0 vs i>0)

	} // END for (i; loop over secondary particles)

	  //	std::cout << "part5: " << part5 << " \n";  /////////////////////////////////////////////////////////////////////
	  //	std::cout << "part6: " << part6 << " \n";  /////////////////////////////////////////////////////////////////////

	  //	std::cout << "C: part1: " << -part1 << " | part2: " << -part2 << " | part3: " << -part3 << " | part4: " << part4 << " | part5: " << part5 << " | part6: " << part6 << " \n";

	return exp(theta*(-part1 - part2 - part3 + part4 + part5 + part6));
}