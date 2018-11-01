
#include "Header.h"

using namespace voro; // kvuli zavedeni trid (napr kvuli pouziti slova pre_container)
					  
// --------------------------------------------------------------------------------------------------------------------------------------------
// vyhledani castice 

// jak urcit souradnice bodu zname-li jeho ID? 
//		- fc positions - vytvoreni vektoru (datove naroky)(transformace dat)
//		- fc findpos - hleda cyklem polohu bodu v pameti (casove naroky)(vyhodnejsi kvuli particle_order class)

/*
// fce positions nacita predem neupresneny pocet bodu do vektoru; body jsou ulozeny v souboru txt ve formatu kdy na
// kazdem radku jsou udaje tykajici se jednoho bodu, a to ID cislo, souradnice x, y, z
void positions(std::vector<float> &v)
{
	// [out]	v		vector for storing the particles id and coordinates.

	int j = 0; // pocita pocet bodu
	int a;
	float b;
	FILE *f;  // deklarace ukazatelu
	f = fopen("data.txt", "r");
	fscanf(f, "%d", &a);

	while (feof(f) == 0) {
		fscanf(f, "%f", &b);
		v.push_back(b);
		fscanf(f, "%f", &b);
		v.push_back(b);
		fscanf(f, "%f", &b);
		v.push_back(b);
		fscanf(f, "%d", &a);
		j++;
	}

	fclose(f);

	// v.size() - vraci velikost vektoru; souradnice vektoru se cisluji od 0!!!
	// v[4k] - id (k+1) bodu, k in (0,1,...)
	// v[4k+1], v[4k+2], v[4k+3] - x, y, z souradnice (k+1)teho bodu 
} 
// NEPOUZITA*/

// fce find_pos najde umisteni castice s danym ID v containeru; varianty pro container a container_poly
void find_pos(int &rijk, int &rq, const int &rid, const voro::container * const pcon)
{
	// [in]		id		id of the searched particle.
	// [in]		con		the container with stored particles.
	// [out]	ijk		the block into which the particle was placed.
	// [out]	q		the position within the block where the particle was placed.

	int i, j;
	bool ch = 0;
	int nboxes = (*pcon).nxyz; // number of computational boxes in container. 
	rijk = 0; rq = 0; // inicializace, aby se vyhlo pripadnym chybam s pameti
	for (i = 0; i < nboxes; i++) {
		for (j = 0; j < (*pcon).co[i]; j++) {
			if ((*pcon).id[i][j] == rid) {
				rijk = i;
				rq = j;
				ch = 1;
				break;
			}
		}
		if (ch = 0) { std::cout << "Particle NOT found! \n"; }
		// computation of coordinates: x = con.p[ijk][3*q]; y = con.p[ijk][3*q+1]; z = con.p[ijk][3*q+2]; 
	}
}

void find_pos(int &rijk, int &rq, const int &rid, const voro::container_poly * const pcon)
{
	// [in]		id		id of the searched particle.
	// [in]		con		the container with stored particles.
	// [out]	ijk		the block into which the particle was placed.
	// [out]	q		the position within the block where the particle was placed.

	int i, j;
	bool ch = 0;
	int nboxes = (*pcon).nxyz; // number of computational boxes in container. 
	rijk = 0; rq = 0; // inicializace, aby se vyhlo pripadnym chybam s pameti
	for (i = 0; i < nboxes; i++) {
		for (j = 0; j < (*pcon).co[i]; j++) {
			if ((*pcon).id[i][j] == rid) {
				rijk = i;
				rq = j;
				ch = 1;
				break;
			}
		}
		if (ch = 0) { std::cout << "Particle NOT found! \n"; }
		// computation of coordinates: x = con.p[ijk][3*q]; y = con.p[ijk][3*q+1]; z = con.p[ijk][3*q+2]; 
	}
}

// najde castici podle jejiho poradoveho cisla, napr. pro 1000 castic je mozne hledat castici s poradovym cislem 1 az 1000; varianty pro container a container_poly
bool find_part(int &ijk, int &q, const int &no, voro::container * const pcon)  // container nemuze byt const kvuli total_particles
{
	// [in]		no		order number of searched particle (order wrt computational boxes).
	// [in]		con		the container with stored particles.
	// [out]	ijk		the block into which the particle was placed.
	// [out]	q		the position within the block where the particle was placed.
	int cit = 0;
	int nu = 0;
	ijk = 0; q = 0; // inicializace

	if (no > (*pcon).total_particles()) { std::cout << "ERROR: order number is higher than total number of particles \n"; return false; }
	else {
		do { nu = nu + (*pcon).co[cit++]; } while (nu < no);
		ijk = cit - 1;
		q = (*pcon).co[cit - 1] - (nu - no + 1);
	}
	return true;
}

bool find_part(int &ijk, int &q, const int &no, voro::container_poly * const pcon)  // container nemuze byt const kvuli total_particles
{
	// [in]		no		order number of searched particle (order wrt computational boxes).
	// [in]		con		the container with stored particles.
	// [out]	ijk		the block into which the particle was placed.
	// [out]	q		the position within the block where the particle was placed.
	int cit = 0;
	int nu = 0;
	ijk = 0; q = 0; // inicializace

	if (no > (*pcon).total_particles()) { std::cout << "ERROR: order number is higher than total number of particles \n"; return false; }
	else {
		do { nu = nu + (*pcon).co[cit++]; } while (nu < no);
		ijk = cit - 1;
		q = (*pcon).co[cit - 1] - (nu - no + 1);
	}
	return true;
}


// --------------------------------------------------------------------------------------------------------------------------------------------
// smazani castice z containeru

// fce erase vymaze danou castici z containeru = z jeho dat struktur id,p,co ; 
// ID smazane castice by fce mela oznacit za znovu pouzitelne - toto ID se ulozi do vektoru fid; varianty pro container a container_poly
bool erase(const int &rijk, int &rq, std::vector<int> *pfid, container *pcon)
{
	// [in,out]		con		the container with stored particles.
	// [in,out]		fid		vector with available id numbers.
	// [in]			ijk,q	the position of the deleted particle.

	if (rijk > ((*pcon).nxyz - 1) || rq > ((*pcon).co[rijk] - 1)) { std::cout << "ERROR: invalid particle"; return false; }
	(*pfid).push_back((*pcon).id[rijk][rq]); // ulozeni ID mazane castice do fid, aby mohlo byt znovu pouzito

	while (rq < ((*pcon).co[rijk] - 1)) {
		(*pcon).id[rijk][rq] = (*pcon).id[rijk][rq + 1];
		(*pcon).p[rijk][3 * rq] = (*pcon).p[rijk][3 * (rq + 1)];
		(*pcon).p[rijk][3 * rq + 1] = (*pcon).p[rijk][3 * (rq + 1) + 1];
		(*pcon).p[rijk][3 * rq + 2] = (*pcon).p[rijk][3 * (rq + 1) + 2];
		rq++;
	}    // v arrays id a p dojde k posunu = vynechani informace o mazane castici
	(*pcon).co[rijk] -= 1;  // zmensi pocet castic v boxu ijk o jedna

							// mem (naalokovana pamet) zustane nezmenena; ... jeste neco chybi ???
	return true;
}  

bool erase(const int &rijk, int &rq, std::vector<int> *pfid, container_poly *pcon)
{
	// [in,out]		con		the container with stored particles.
	// [in,out]		fid		vector with available id numbers.
	// [in]			ijk,q	the position of the deleted particle.

	if (rijk > ((*pcon).nxyz - 1) || rq > ((*pcon).co[rijk] - 1)) { std::cout << "ERROR: invalid particle"; return false; }
	(*pfid).push_back((*pcon).id[rijk][rq]); // ulozeni ID mazane castice do fid, aby mohlo byt znovu pouzito

	while (rq < ((*pcon).co[rijk] - 1)) {
		(*pcon).id[rijk][rq] = (*pcon).id[rijk][rq + 1];
		(*pcon).p[rijk][4 * rq] = (*pcon).p[rijk][4 * (rq + 1)];
		(*pcon).p[rijk][4 * rq + 1] = (*pcon).p[rijk][4 * (rq + 1) + 1];
		(*pcon).p[rijk][4 * rq + 2] = (*pcon).p[rijk][4 * (rq + 1) + 2];
		(*pcon).p[rijk][4 * rq + 3] = (*pcon).p[rijk][4 * (rq + 1) + 3];
		rq++;
	}    // v arrays id a p dojde k posunu = vynechani informace o mazane castici
	(*pcon).co[rijk] -= 1;  // zmensi pocet castic v boxu ijk o jedna

							// mem (naalokovana pamet) zustane nezmenena; ... jeste neco chybi ???
	return true;
}

// druha varianta fce erase, tentokrat bez ukladani do vektoru fid; varianty pro container a container_poly
bool erase(const int &rijk, int &rq, container *pcon)
{
	// [in,out]		con		the container with stored particles.
	// [in]			ijk,q	the position of the deleted particle.

	if (rijk > ((*pcon).nxyz - 1) || rq > ((*pcon).co[rijk] - 1)) { std::cout << "ERROR: invalid particle"; return false; }

	while (rq < ((*pcon).co[rijk] - 1)) {
		(*pcon).id[rijk][rq] = (*pcon).id[rijk][rq + 1];
		(*pcon).p[rijk][3 * rq] = (*pcon).p[rijk][3 * (rq + 1)];
		(*pcon).p[rijk][3 * rq + 1] = (*pcon).p[rijk][3 * (rq + 1) + 1];
		(*pcon).p[rijk][3 * rq + 2] = (*pcon).p[rijk][3 * (rq + 1) + 2];
		rq++;
	}    // v arrays id a p dojde k posunu = vynechani informace o mazane castici
	(*pcon).co[rijk] -= 1;  // zmensi pocet castic v boxu ijk o jedna

							// mem (naalokovana pamet) zustane nezmenena; ... jeste neco chybi ???
	return true;
}

bool erase(const int &rijk, int &rq, container_poly *pcon)
{
	// [in,out]		con		the container with stored particles.
	// [in]			ijk,q	the position of the deleted particle.

	if (rijk > ((*pcon).nxyz - 1) || rq > ((*pcon).co[rijk] - 1)) { 
		//std::cout << "ERROR: invalid particle"; 
		std::cout << "ERROR: invalid particle " << rijk << " " << rq << "\n";
		return false; }

	while (rq < ((*pcon).co[rijk] - 1)) {
		(*pcon).id[rijk][rq] = (*pcon).id[rijk][rq + 1];
		(*pcon).p[rijk][4 * rq] = (*pcon).p[rijk][4 * (rq + 1)];
		(*pcon).p[rijk][4 * rq + 1] = (*pcon).p[rijk][4 * (rq + 1) + 1];
		(*pcon).p[rijk][4 * rq + 2] = (*pcon).p[rijk][4 * (rq + 1) + 2];
		(*pcon).p[rijk][4 * rq + 3] = (*pcon).p[rijk][4 * (rq + 1) + 3];
		rq++;
	}    // v arrays id a p dojde k posunu = vynechani informace o mazane castici
	(*pcon).co[rijk] -= 1;  // zmensi pocet castic v boxu ijk o jedna

							// mem (naalokovana pamet) zustane nezmenena; ... jeste neco chybi ???
	return true;
}

// --------------------------------------------------------------------------------------------------------------------------------------------
// determining neighbours

// fce are_neighbors pro danou bunku c zjisti, zda je dana castice jejim sousedem; varianty pro container a container_poly
bool are_neighbors(voronoicell_neighbor &rc, const int &rijk, const int &rq, const voro::container * const pcon)
{
	// [in]		c		the first considered cell.
	// [in]		ijk		the block of particle of the second considered cell.
	// [in]		q		the position within the block of particle of the second considered cell.
	// [in]		con		the container with stored particles.

	std::vector<int> sous;  // vektor pro prirazeni ID sousedu
	unsigned int i;
	rc.neighbors(sous);		// spocte ID sousedu do tohoto vektoru
	for (i = 0; i < sous.size(); i++) { // loop over the neighbors 
		if (sous[i] == (*pcon).id[rijk][rq]) { return true; }
	}
	return false;
}

bool are_neighbors(voronoicell_neighbor &rc, const int &rijk, const int &rq, const voro::container_poly * const pcon)
{
	// [in]		c		the first considered cell.
	// [in]		ijk		the block of particle of the second considered cell.
	// [in]		q		the position within the block of particle of the second considered cell.
	// [in]		con		the container with stored particles.

	std::vector<int> sous;  // vektor pro prirazeni ID sousedu
	unsigned int i;
	rc.neighbors(sous);		// spocte ID sousedu do tohoto vektoru
	for (i = 0; i < sous.size(); i++) { // loop over the neighbors 
		if (sous[i] == (*pcon).id[rijk][rq]) { return true; }
	}
	return false;
}

// fce secondary vezme souseda sekundarni castice a pokud se jedna taktez o sekundarni castici vrati jeho polohu ijk,q a true - NEPOUZITA
bool secondary(const int cid, const std::vector<int> sec, const voro::container * const pcon, int &ijk, int &q)
{
	// [in]		cid		id of the current particle under consideration.
	// [in]		sec		the vector of secondary particles (neighbors of primary particle).
	// [in]		con		the container with stored particles.
	// [out]	ijk,q	position of the neighbor.

	unsigned int i;
	for (i = 0; i < sec.size() / 2; i++) {
		if (cid == (*pcon).id[sec[2 * i]][sec[2 * i + 1]]) { ijk = sec[2 * i]; q = sec[2 * i + 1]; return true; }
	}
	return false;
}

// fce terciary vezme souseda sekundarni castice, a pokud se neshoduje ani s primarni ani s zadnou sekundarni castici oznaci jej jako terciarni castici;  varianty pro container a container_poly
bool terciary(const int cid, const int pid, const std::vector<int> sec, const voro::container * const pcon)
{
	// [in]		cid		id of the current particle under consideration.
	// [in]		pid		id of the primary particle (added particle).
	// [in]		sec		the vector of secondary particles (neighbors of primary particle).
	// [in]		con		the container with stored particles.

	unsigned int i;
	if (cid == pid) { return false; }
	for (i = 0; i < sec.size() / 2; i++) {
		if (cid == (*pcon).id[sec[2 * i]][sec[2 * i + 1]]) { return false; }
	}
	return true;
}

bool terciary(const int cid, const int pid, const std::vector<int> sec, const voro::container_poly * const pcon)
{
	// [in]		cid		id of the current particle under consideration.
	// [in]		pid		id of the primary particle (added particle).
	// [in]		sec		the vector of secondary particles (neighbors of primary particle).
	// [in]		con		the container with stored particles.

	unsigned int i;
	if (cid == pid) { return false; }
	for (i = 0; i < sec.size() / 2; i++) {
		if (cid == (*pcon).id[sec[2 * i]][sec[2 * i + 1]]) { return false; }
	}
	return true;
}

// when the primary particle was deleted we use this case:
bool terciary(const int cid, const std::vector<int> sec, const voro::container * const pcon)
{
	// [in]		cid		id of the current particle under consideration.
	// [in]		sec		the vector of secondary particles (neighbors of primary particle).
	// [in]		con		the container with stored particles.

	unsigned int i;
	for (i = 0; i < sec.size() / 2; i++) {
		if (cid == (*pcon).id[sec[2 * i]][sec[2 * i + 1]]) { return false; }
	}
	return true;
}

bool terciary(const int cid, const std::vector<int> sec, const voro::container_poly * const pcon)
{
	// [in]		cid		id of the current particle under consideration.
	// [in]		sec		the vector of secondary particles (neighbors of primary particle).
	// [in]		con		the container with stored particles.

	unsigned int i;
	for (i = 0; i < sec.size() / 2; i++) {
		if (cid == (*pcon).id[sec[2 * i]][sec[2 * i + 1]]) { return false; }
	}
	return true;
}


// -----------------------------------------------------------------------------------------------------------------------------------
// uprava vektoru

// fce identical porovna dva vektory (slozky kazdeho vektoru jsou navzajem ruzne) a rozhodne zda se lisi ci ne - NEPOUZITA
bool identical(const std::vector<int> &ra, const std::vector<int> &rb)
{
	// [in]		a,b		two vectors to be compared.

	unsigned int i, j;
	bool shoda;
	for (i = 0; i < ra.size(); i++) {
		shoda = 0;
		for (j = 0; j < rb.size(); j++) {
			if (ra[i] == rb[j]) { shoda = 1; }
		}
		if (shoda == 0) { return false; }
	}
	return true;
}

// function merge merges vectors containing information about secondary particles
void merge(std::vector<int> &rna, std::vector<int> &rnb, std::vector<int> &ra, std::vector<int> &rb, int k,
	std::vector<int> &raa, std::vector<int> &rbb, std::vector<double> &raaa, std::vector<double> &rbbb)
{
	// [in]			na		the first vector containing neighborhood information.
	// [in]			nb		the second vector containing neighborhood information.
	// [in,out]		a		the first vector of secondary particles, the second one is merged into this vector.
	// [in]			b		the second vector of secondary particles.
	// [in]			k		the number of secondary particles which are "out" in the first vector.
	// [in,out]		aa		the vector containing "in"/"out" information for the the first vector of secondary particles (used for merging).
	// [in]			bb		the vector containing "in"/"out" information for the the second vector of secondary particles.
	// [in,out]		aaa		the vector containing the real coordinates of "out" particles from the first vector (used for merging).
	// [in]			bbb		the vector containing the real coordinates of "out" particles from the second vector.

	// example of usage: merge(neigh_add, neigh_del, sr_add, sr_del, k1_add, sio_add, sio_del, sap_add, sap_del);

	unsigned int i, j;
	unsigned int length = rna.size();
	bool shoda;

	for (i = 0; i < rnb.size(); i++) {  // go through the second vector
		shoda = 0;
		for (j = 0; j < length; j++) {  // try if actually considered particle of the second vector is inside the first vector
			if (rnb[i] == rna[j]) { shoda = 1; break; }
		}
		if (shoda == 0) {				// if not change vectors
			rna.push_back(rnb[i]);
			ra.push_back(rb[2 * i]);	ra.push_back(rb[2 * i + 1]);									// merge vectors sr
			if (rbb[i] > 0) {																		// merge vectors sio
				k++;
				raa.push_back(k);
				// if value in sio is positive then merge vectors sap
				raaa.push_back(rbbb[3 * (rbb[i] - 1)]); raaa.push_back(rbbb[3 * (rbb[i] - 1) + 1]); raaa.push_back(rbbb[3 * (rbb[i] - 1) + 2]); // merge vectors sap
			}
			else { raa.push_back(0); }

		}
	}
}


// merge of two integer vectors
void merge(std::vector<int> &rna, std::vector<int> &rnb)
{
	// [in,out]		na		the first vector 
	// [in]			nb		the second vector 

	unsigned int i, j;
	unsigned int length = rna.size();
	bool shoda;

	for (i = 0; i < rnb.size(); i++) {  // go through the second vector
		shoda = 0;
		for (j = 0; j < length; j++) {  // try if actually considered particle of the second vector is inside the first vector
			if (rnb[i] == rna[j]) { shoda = 1; break; }
		}
		if (shoda == 0) {				// if not change vectors
			rna.push_back(rnb[i]);
		}
	}
}


// merge vector and integer
bool merge(std::vector<int> &rna, int &rnb)
{
	// [in,out]		na		the first vector 
	// [in]			nb		the integer 

	unsigned int i, j;
	unsigned int length = rna.size();
	bool shoda = 0;

	for (j = 0; j < length; j++) {  // try if the value nb is inside the vector
		if (rnb == rna[j]) { shoda = 1; break; }
	}
	if (shoda == 0) {				// if not add it to the end of the vector
		rna.push_back(rnb);
		return true;
	}
	return false;
}

// fcs vector_dif give us elements which are not common in the given two vectors
//	two functions = two different ways how to do it
void vector_dif(std::vector<int> &rna, std::vector<int> &rnb)
{
	// [in,out]		na		the first vector 
	// [in]			nb		the second vector 

	int length = rna.size()-1;
	int i,j;

	for (i = length; i >= 0; i--) {
		for (j = (rnb.size() - 1); j >= 0; j--) {
			if (rna[i] == rnb[j]) {				// the pair was found
				rna.erase(rna.begin() + i);		// erase the element in the first vector
				rnb.erase(rnb.begin() + j);		// erase the element in the second vector
				break;							// end the second for loop
			}
		}
	}
}

// returns vector of elements of na which are not common to nb
void vector_dif(std::vector<int> &rna, std::vector<int> &rnb, std::vector<int> &rnc)
{
	// [in]			na		the first vector 
	// [in]			nb		the second vector 

	//int length = rna.size() - 1;
	int i, j;
	bool shoda = false;

	//rnc.clear();

	for (i = 0; i < rna.size(); i++) {
		shoda = false;
		for (j = 0; j < rnb.size(); j++) {
			if (rna[i] == rnb[j]) { shoda = true; break; }			// the pair was found	
		}
		if (shoda == false) { rnc.push_back(rna[i]); }
	}

	//for (i = 0; i<rnb.size(); i++) {
	//	shoda = false;
	//	for (j = 0; j < rna.size(); j++) {
	//		if (rna[i] == rnb[j]) { shoda = true; }			// the pair was found	
	//	}
	//	if (shoda == false) { rnc.push_back(rnb[i]); }
	//}
}


void sec_neigh(std::vector<int> &rcells, std::vector<int> &rsec)
{
	//	[in]	con		container with stored particles
	//	[in]	rcells	group of cells to which we want to determine secondary neighborhood (neighbours of these cells not included in the group)
	//	[out]	rsec	vector of secondary cells

	voronoicell_neighbor c;

}

// ---------------------------------------------------------------------------------------------------------------------------------------------
// hledani minima a maxima

// fce min_max usporada dve hodnoty a,b tak, aby a >= b
void min_max(double &ra, double &rb)
{
	// [in,out]		a		the first value to be sorted, after sorting a >= b.
	// [in,out]		b		the second value to be sorted.

	double c = 0;
	if (rb > ra) { c = ra; ra = rb; rb = c; } // prohozeni hodnot
}




// fce min_max usporada dve hodnoty a,b tak, aby a >= b
void min_max(int &ra, int &rb)
{
	// [in,out]		a		the first value to be sorted, after sorting a >= b.
	// [in,out]		b		the second value to be sorted.

	int c = 0;
	if (rb > ra) { c = ra; ra = rb; rb = c; } // prohozeni hodnot
}

// fce min_max usporada tri hodnoty a,b,c tak, aby c >= b >= a
void min_max(double &ra, double &rb, double &rc)
{
	// [in]		a,b,c		three values to be ordered as a <= b <= c

	double min, max;

	if (ra < rb) { min = ra; max = rb; }
	else { min = rb; max = ra; }
	if (rc < min) { ra = rc; rb = min; rc = max; return; }
	if (rc > max) { ra = min; rb = max; return; }

	ra = min; rb = rc; rc = max;
}

// fce min_max usporada tri hodnoty a,b,c tak, aby c >= b >= a
void min_max(int &ra, int &rb, int &rc)
{
	// [in]		a,b,c		three values to be ordered as a <= b <= c

	int min, max;

	if (ra < rb) { min = ra; max = rb; }
	else { min = rb; max = ra; }
	if (rc < min) { ra = rc; rb = min; rc = max; return; }
	if (rc > max) { ra = min; rb = max; return; }

	ra = min; rb = rc; rc = max;
}


double abs_val(double val) { 
	if (val > 0) { return val; } 
	else { return (val - 2 * val); } 
}

// ------------------------------------------------------------------------------------------------------------------------------------------
// teziste, normaly, vzdalenosti sten a bodu

// fce barycentrum spocte barycentrum dvou bunek
bool barycentrum(double &rx, double &ry, double &rz, double &rnx, double &rny, double &rnz, double ra, double rb)
{
	// [in]		x,y,z		the coordinates of the centroid of the first particle.
	// [in]		nx,ny,nz	the coordinates of the centroid of the second particle.
	// [in]		a			the volume of the first cell.
	// [in]		b			the volume of the second particle.

	double tx, ty, tz; // coordinates of barycenter
	double dec = 10000000000000; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI

	tx = (ra*rx + rb*rnx) / (ra + rb);  // tx = (a*x + b*(nx+sx))/(a+b)
	tx = round(tx*dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	if (tx > 1 || tx < 0) { return false; } // teziste je mimo okno --> nepocitej fci V2
	ty = (ra*ry + rb*rny) / (ra + rb);  // ty = (a*y + b*(ny+sy))/(a+b)
	ty = round(ty*dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	if (ty > 1 || ty < 0) { return false; } // teziste je mimo okno --> nepocitej fci V2
	tz = (ra*rz + rb*rnz) / (ra + rb);  // tz = (a*z + b*(nz+sz))/(a+b)
	tz = round(tz*dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	if (tz > 1 || tz < 0) { return false; } // teziste je mimo okno --> nepocitej fci V2

	return true;  // barycenter is inside the window
}

// fce bar_coor urci souradnice barycentra - NEPOUZITA, zbytecna
void bar_coor(double &rx, double &ry, double &rz, double &rnx, double &rny, double &rnz, double &ra, double &rb) // NEPOUZITA
{
	// [in, out]	x,y,z		the coordinates of the first particle.
	// [in]			nx,ny,nz	the coordinates of the second particle.
	// [in]			a			the volume of the first cell.
	// [in]			b			the volume of the second particle.

	// zkusme to bez zaokrouhleni
	
	rx = (ra*rx + rb*rnx) / (ra + rb);  // tx = (a*x + b*(nx+sx))/(a+b)
	ry = (ra*ry + rb*rny) / (ra + rb);  // ty = (a*y + b*(ny+sy))/(a+b)
	rz = (ra*rz + rb*rnz) / (ra + rb);  // tz = (a*z + b*(nz+sz))/(a+b)
	
}

void face_dist(const unsigned int &rfing, const std::vector<int> &fvert, const double &rx, const double &ry, const double &rz,
	double &rnx, double &rny, double &rnz, voronoicell_neighbor &rc)
{
	// [in]		fing		vector index indicating where is information about given face starting.
	// [in]		fvert		vector containing a list of vertices for each face
	// [in]		x,y,z		the coordinates of the particle.
	//				souradnice musi byt v okne, spocte se normalovy vektor, skutecna poloha souradnic se pouzije az potom - pricte se k ni 2* tento vektor
	// [out]	nx,ny,nz	the normal vector of the face with the length equal to the face distance
	// [in]		c			voronoi cell.

	double a1, a2, a3, b1, b2, b3, c1, c2, c3, u1, u2, u3, v1, v2, v3, a, b, c, d;

	//double dec = 10000000000000; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI

	if (fvert[rfing] < 3) { std::cout << "face_dist: Not enough vertices! \n"; return; }
	// determine 3 points of the given face, compute directions vector and its vector product
	// vertices: A = fvert[rfing+1], B = fvert[rfing+2], C = fvert[rfing+3]
	// coordinates of the 1st: A = [ pts[3*fvert[rfing+1]], pts[3*fvert[rfing+1] + 1], pts[3*fvert[rfing+1] + 2] ]
	a1 = rx + rc.pts[3 * fvert[rfing + 1]] * 0.5; a2 = ry + rc.pts[3 * fvert[rfing + 1] + 1] * 0.5; a3 = rz + rc.pts[3 * fvert[rfing + 1] + 2] * 0.5;
	//std::cout << a1 << " " << a2 << " " << a3 << " " << '\n';  //////////////////////////////////////////////////////////
	b1 = rx + rc.pts[3 * fvert[rfing + 2]] * 0.5; b2 = ry + rc.pts[3 * fvert[rfing + 2] + 1] * 0.5; b3 = rz + rc.pts[3 * fvert[rfing + 2] + 2] * 0.5;
	//std::cout << b1 << " " << b2 << " " << b3 << " " << '\n';  //////////////////////////////////////////////////////////
	c1 = rx + rc.pts[3 * fvert[rfing + 3]] * 0.5; c2 = ry + rc.pts[3 * fvert[rfing + 3] + 1] * 0.5; c3 = rz + rc.pts[3 * fvert[rfing + 3] + 2] * 0.5;
	//std::cout << c1 << " " << c2 << " " << c3 << " " << '\n';  //////////////////////////////////////////////////////////
	// directions vectors: u = B-A , v = C-A
	u1 = b1 - a1; u2 = b2 - a2; u3 = b3 - a3;
	v1 = c1 - a1; v2 = c2 - a2; v3 = c3 - a3;
	// vector product: w = u*v  --> coefficients a,b,c in the plane equation
	// a = u2*v3-v2*u3 , b = u3*v1-v3*u1 , c = u1*v2-v1*u2
	a = u2*v3 - v2*u3; b = u3*v1 - v3*u1; c = u1*v2 - v1*u2;
	// computing coefficient d by plugging in the vertex coordinates
	d = a*c1 + b*c2 + c*c3;
	// we arrive to the plane equation ax+by+cz=d

	// zaokrouhleni koeficientu rovnice primky
	//a = round(a * dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	//b = round(b * dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	//c = round(c * dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	//d = round(d * dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	
	// translation: [rx,ry,rz] will be new origin - translated plane a(x+rx)+b(y+ry)+c(z+rz)=d
	// chenge of coordinated results in the plane equation ax+by+cz=d-(arx+bry+crz), i.e. ax+by+cz=d*, where d*=d-(arx+bry+crz) ... b1=-d*
	b1 = a*rx + b*ry + c*rz - d;
	// distance: |a*rx+b*ry+c*rx+d|/sqrt(a^2+b^2+c^2)
	//b2 = a*rx + b*ry + c*rz + d;
	b3 = pow(a, 2) + pow(b, 2) + pow(c, 2);

	//std::cout << a << " " << b << " " << c << " " << d << '\n'; /////////////////////////////////////////////////////////
	//std::cout << b1 << " " << b3 << '\n';  //////////////////////////////////////////////////////////////////////////////
	//std::cout << rx - a*(b1 / b3) << " " << ry - b*(b1 / b3) << " " << rz - c*(b1 / b3) << " " << '\n'; /////////////////

	// rnx = (rx - a*b1 / b3) - rx; rny = (ry - b*b1 / b3) - ry; rnz = (rz - b*b1 / b3) - rz; -->
	rnx = -(a*b1 / b3); rny = -(b*b1 / b3); rnz = -(c*b1 / b3); // not unit normal vector

	//rnx = round(rnx * dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	//rny = round(rny * dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	//rnz = round(rnz * dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
								  //if (b2 < 0) { b2 = -b2; }

								  //return b2 / sqrt(b3);

}

void face_dist(const unsigned int &rfing, const std::vector<int> &fvert, const double &rx, const double &ry, const double &rz,
	const double &tx, const double &ty, const double &tz, double &rnx, double &rny, double &rnz, voronoicell_neighbor &rc)
{
	// [in]		fing		vector index indicating where is information about given face starting.
	// [in]		fvert		vector containing a list of vertices for each face
	// [in]		rx,ry,rz	the coordinates of the generator.
	//				souradnice musi byt v okne, spocte se normalovy vektor, skutecna poloha souradnic se pouzije az potom - pricte se k ni 2* tento vektor
	// [in]		tx,ty,tz	the coordinates of the barycenter
	// [out]	nx,ny,nz	the normal vector of the face with the length equal to the face distance
	// [in]		c			voronoi cell.

	double a1, a2, a3, b1, b2, b3, c1, c2, c3, u1, u2, u3, v1, v2, v3, a, b, c, d;

	//double dec = 10000000000000; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI

	if (fvert[rfing] < 3) { std::cout << "face_dist: Not enough vertices! \n"; return; }
	// determine 3 points of the given face, compute directions vector and its vector product
	// vertices: A = fvert[rfing+1], B = fvert[rfing+2], C = fvert[rfing+3]
	// coordinates of the 1st: A = [ pts[3*fvert[rfing+1]], pts[3*fvert[rfing+1] + 1], pts[3*fvert[rfing+1] + 2] ]
	a1 = rx + rc.pts[3 * fvert[rfing + 1]] * 0.5; a2 = ry + rc.pts[3 * fvert[rfing + 1] + 1] * 0.5; a3 = rz + rc.pts[3 * fvert[rfing + 1] + 2] * 0.5;
	//std::cout << a1 << " " << a2 << " " << a3 << " " << '\n';  //////////////////////////////////////////////////////////
	b1 = rx + rc.pts[3 * fvert[rfing + 2]] * 0.5; b2 = ry + rc.pts[3 * fvert[rfing + 2] + 1] * 0.5; b3 = rz + rc.pts[3 * fvert[rfing + 2] + 2] * 0.5;
	//std::cout << b1 << " " << b2 << " " << b3 << " " << '\n';  //////////////////////////////////////////////////////////
	c1 = rx + rc.pts[3 * fvert[rfing + 3]] * 0.5; c2 = ry + rc.pts[3 * fvert[rfing + 3] + 1] * 0.5; c3 = rz + rc.pts[3 * fvert[rfing + 3] + 2] * 0.5;
	//std::cout << c1 << " " << c2 << " " << c3 << " " << '\n';  //////////////////////////////////////////////////////////
	// directions vectors: u = B-A , v = C-A
	u1 = b1 - a1; u2 = b2 - a2; u3 = b3 - a3;
	v1 = c1 - a1; v2 = c2 - a2; v3 = c3 - a3;
	// vector product: w = u*v  --> coefficients a,b,c in the plane equation
	// a = u2*v3-v2*u3 , b = u3*v1-v3*u1 , c = u1*v2 - v1*u2
	a = u2*v3 - v2*u3; b = u3*v1 - v3*u1; c = u1*v2 - v1*u2;
	// computing coefficient d by plugging in the vertex coordinates
	d = a*c1 + b*c2 + c*c3;

	// zaokrouhleni koeficientu rovnice primky
	//a = round(a * dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	//b = round(b * dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	//c = round(c * dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI
	//d = round(d * dec) / dec; //ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI ZAOKROUHLENI

	b1 = a*tx + b*ty + c*tz - d;
	// distance: |a*rx+b*ry+c*rx+d|/sqrt(a^2+b^2+c^2)
	//b2 = a*rx + b*ry + c*rz + d;
	b3 = pow(a, 2) + pow(b, 2) + pow(c, 2);

	//double dist = abs(b1) / sqrt(b3);
	//if (dist > 1) {
	//	std::cout << dist << " ; " << rx << " " << ry << " " << rz << " ; " << tx << " " << ty << " " << tz << '\n';
	//}

	//std::cout << a << " " << b << " " << c << " " << d << '\n'; /////////////////////////////////////////////////////////
	//std::cout << b1 << " " << b3 << '\n';  //////////////////////////////////////////////////////////////////////////////
	//std::cout << rx - a*(b1 / b3) << " " << ry - b*(b1 / b3) << " " << rz - c*(b1 / b3) << " " << '\n'; /////////////////

	// rnx = (rx - a*b1 / b3) - rx; rny = (ry - b*b1 / b3) - ry; rnz = (rz - b*b1 / b3) - rz; -->
	rnx = -(a*b1 / b3); rny = -(b*b1 / b3); rnz = -(c*b1 / b3); // not unit normal vector

										
}


// fc real_coo for given pair of generators (aasumed to be neighbors) returns true values of coordinates of the second generator
// true = means that coordinates can be outside the window; function is suitable only for beta < 1/4 (beta = hardcore parameter)
// warning: using beta > 1/4 is not recomended, can be unstable; therefore this fc is sufficient and can be used instead of fc face_dist (in Voronoi case) and for Laguerre too
void real_coo(double x, double y, double z, double &xx, double &yy, double &zz)
{
	// [in]			x,y,z		the coordinates of the first particle.
	// [in/out]		xx,yy,zz	the coordinates of the second particle (that in the window on input, the real ones on the output).

	// x coordinate:
	if (abs(x - xx) < 1 / 2) {} // vzd bodu mensi nez 1/2 --> xx je v okne
	else {
		if (x - xx > 1 / 2) { xx = xx + 1; }
		if (x - xx < -1 / 2) { xx = xx - 1; }
	}
	// y coordinate:
	if (abs(y - yy) < 1 / 2) {} // vzd bodu mensi nez 1/2 --> yy je v okne
	else {
		if (y - yy > 1 / 2) { yy = yy + 1; }
		if (y - yy < -1 / 2) { yy = yy - 1; }
	}
	// z coordinate:
	if (abs(z - zz) < 1 / 2) {} // vzd bodu mensi nez 1/2 --> zz je v okne
	else {
		if (z - zz > 1 / 2) { zz = zz + 1; }
		if (z - zz < -1 / 2) { zz = zz - 1; }
	}
}

double point_dist(double &rx, double &ry, double &rz, double &rnx, double &rny, double &rnz)
{
	// [in]		x,y,z		the coordinates of the first particle.
	// [in]		nx,ny,nz	the coordinates of the second particle.

	double ni, nj, nk;
	// vzdalenost dvou bodu: sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
	ni = pow(rx - rnx, 2);      // double pow (double base, double exponent);
	nj = pow(ry - rny, 2);
	nk = pow(rz - rnz, 2);

	return sqrt(ni + nj + nk);
}

double h_maximum(unsigned int &n, std::vector<int> &vert, voronoicell_neighbor &d, double &x, double &y, double &z)
{
	// [in]		n			number of neighbors
	// [in]		vert		vertices ordered by faces
	// [in]		d			the cell under consideration
	// [in]		x,y,z		its coordinates
	// [out]	h_max		the maximum distance from centre to face

	unsigned int j, fng;
	double xn, yn, zn;
	double dist, h_max;

	h_max = 0;
	fng = 0;
	for (j = 0; j < n; j++) {
		face_dist(fng, vert, x, y, z, xn, yn, zn, d);		// {PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
		xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;	// {PER} urcim skutecne souradnice tohoto souseda za stenou

		dist = point_dist(x, y, z, xn, yn, zn)/2;		// tj dist se nyni rovna norme vektoru (xn,yn,zn) ktery vyjde primo z fce face_dist
		if (dist > h_max) { h_max = dist; }

		fng = fng + vert[fng] + 1;
	}

	return h_max;
}

double h_minimum(unsigned int &n, std::vector<int> &vert, voronoicell_neighbor &d, double &x, double &y, double &z)
{
	// [in]		n			number of neighbors
	// [in]		vert		vertices ordered by faces
	// [in]		d			the cell under consideration
	// [in]		x,y,z		its coordinates
	// [out]	h_min		the maximum distance from centre to face

	unsigned int j, fng;
	double xn, yn, zn;
	double dist, h_min;

	h_min = 2;
	fng = 0;
	for (j = 0; j < n; j++) {
		face_dist(fng, vert, x, y, z, xn, yn, zn, d);		// {PER} urcim vzdalenost steny kterou bunka sdili s aktualnim sousedem		
		xn = x + 2 * xn; yn = y + 2 * yn; zn = z + 2 * zn;	// {PER} urcim skutecne souradnice tohoto souseda za stenou

		dist = point_dist(x, y, z, xn, yn, zn) / 2;		// tj dist se nyni rovna norme vektoru (xn,yn,zn) ktery vyjde primo z fce face_dist
		if (dist < h_min) { h_min = dist; }

		fng = fng + vert[fng] + 1;
	}

	return h_min;
}

void h_fcs(voronoicell_neighbor &d, double &x, double &y, double &z, double &xb, double &yb, double &zb, double &h_max, double &h_min)
{
	// [in]		d				the cell under consideration
	// [in]		x,y,z			the coordinates of its generator (enable to compute the equation of plane of face)
	// [in]		xb,yb,zb		the coordinates of its generator/barycenter
	// [out]	h_max, h_min	the maximal and minimal distance from centre to face

	int j;
	unsigned int fng;
	double xn, yn, zn;
	double dist;
	std::vector<int> vert;

	//d.neighbors(neigh);
	d.face_vertices(vert); 

	h_max = 0;
	h_min = 200000;
	fng = 0;
	for (j = 0; j < d.number_of_faces(); j++) {
		face_dist(fng, vert, x, y, z, xb, yb, zb, xn, yn, zn, d);		// urcim normalovy vektor steny od generatoru/teziste bunky		

		dist = sqrt(pow(xn, 2) + pow(yn, 2) + pow(zn, 2));	// norma tohoto vektoru
	
		if (dist < h_min) { h_min = dist; }
		if (dist > h_max) { h_max = dist; }

		fng = fng + vert[fng] + 1;
	}
}

void volume_min_max(voro::container &rcon)
{
	int i, j;
	double vol, vol_min, vol_max;
	voronoicell c;  // bunka
	
	vol_min = 1;
	vol_max = 0;

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			rcon.compute_cell(c, j, i);
			vol = c.volume();

			if (vol > vol_max) { vol_max = vol; }
			if (vol < vol_min) { vol_min = vol; }

		}
	}
	std::cout << "Max volume: " << vol_max << " , min volume: " << vol_min << "\n";
}

void area_min_max(voro::container &rcon)
{
	int i, j;
	double vol, vol_min, vol_max;
	voronoicell c;  // bunka

	vol_min = 100;
	vol_max = 0;

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			rcon.compute_cell(c, j, i);
			vol = c.surface_area();

			if (vol > vol_max) { vol_max = vol; }
			if (vol < vol_min) { vol_min = vol; }

		}
	}
	std::cout << "Max area: " << vol_max << " , min area: " << vol_min << "\n";

}

// fc find_index finds index of the given value in the given vector (where every element is only once)
// returns -1 if the value was not found
int find_index(std::vector<int> rv, int &ri)
{
	// [in]		v		the vector of which the index is to be found
	// [in]		i		the value we are searching for

	int j;
	j = 0;
	int k;
	k = rv[0];
	int s = rv.size();
	 
	while (k != ri) {
		j++; if ((j+1) > s ) { /*std::cout << "ERROR: find_index - no such value in the given vector! \n";*/ return -1; }
		k = rv[j]; 
	}
	return j;
}



// fc common_edge najde spolecnou hranu dvou sten teze bunky a vrati jeji poradove cislo
// poradove cislo je poradi ve strukture ed, uvazujeme-li pouze hrany, kdy prvni vrchol ma mensi poradove cislo nez druhy vrchol
// vrati-li fce -1, pak steny nemaji spolecnou hranu
int common_edge(voronoicell_neighbor &rc, int &f1, int &f2)
{
	// [in]		c		the voronoi cell with the neighbor information
	// [in]		f1,f2	numbers of the two faces of the cell c

	if (f1 == -1) { return -1; } // unvalid number, e.g. find_index returning -1 means no index found
	if (f2 == -1) { return -1; }

	// note: two faces has no common edge iff appropriate neighboring cells are not neighbors (assumption)
	int e, k, l, p1, p2, v1, v2;
	std::vector<int> ord, vert;
	v1 = -1; v2 = -1;

	rc.face_orders(ord);
	rc.face_vertices(vert);

	p1 = 0;
	for (k = 0; k < f1; k++) {
		p1 = p1 + ord[k];
	}
	p1 = p1 + f1 + 1;			// p1 - prvni vrchol steny f1

	p2 = 0;
	for (k = 0; k < f2; k++) {
		p2 = p2 + ord[k];
	}
	p2 = p2 + f2 + 1;			// p2 - prvni vrchol steny f2

	e = 0;
	for (k = 0; k < ord[f1]; k++) {				//najdi dva spolecne vrcholy techto sten
		for (l = 0; l < ord[f2]; l++) {
			if (vert[p1 + k] == vert[p2 + l]) {
				if (e == 0){ v1 = vert[p1 + k]; e++; }
				else { v2 = vert[p1 + k]; }
			}
		}
	}

	if (e == 0) { return -1; } // vrat -1 pokud nebyl nalezen spolecny vrchol (pozn. pokud steny maji jeden spolecny vrchol, pak musi mit i druhy spolecny vrchol)
	if (v2 == -1) { std::cout << "ERROR: common_edge - druhy vrchol nenalezen \n"; return -1;} // kontrola

	min_max(v1, v2);  // nyni je v1 >= v2

	e = 0;										// najdi poradi teto hrany
	for (k = 0; k < v2; k++) {			// loop pres hrany
		for (l = 0; l < rc.nu[k]; l++) {
			if (k < rc.ed[k][l]) { e++; }
		}
		//e = e + rc.nu[k];
	}
	l = 0;
	e++;
	while (rc.ed[v2][l] != v1) {
		if (v2 < rc.ed[v2][l]) { e++; } 
		l++; 
		if (l >= rc.nu[v2]) { break; std::cout << " ERROR: common_edge - edge not found! \n"; }
	}

	// e je nyni pocet hran, kdy k < ed[k][l] az do hrany [v2,v1] vcetne
	// my ale chceme poradi, nikoliv pocet; poradi bude e-1  - slo by vyrusit e-1 s e++
	// poradi je tedy pocet hran, k < ed[k][l] az do hrany [v2,v1], ale tuto hranu jiz nepocitame
	return e-1;
}


// ----------------------------------------------------------------------------------------------------------------------------------------

// fc point_density returns a vector of counts of points in every set of regular lattice; container and container_poly variants
void point_density(std::vector<int> counts, voro::container &rcon, int &gsi)
{
	// [out]	counts			vector of numbers of points in a given areas of lattice (of the same length as a vector of residuals)
	// [in]		con				the container with stored particles
	// [in]		gsi				the grid size (number of cubes in each direction)

	int i, j;
	int x, y, z;

	counts.resize(gsi*gsi*gsi);
	// loop over all generators -- mod its coordinates -- obtain j,k,l values -- inc appropriate value in the vector
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box


				x = (int)(rcon.p[j][3*i] * gsi);			// coordinates --> i,j,k values
				y = (int)(rcon.p[j][3*i + 1] * gsi);		// (int) : 0.897*10 -> 8 , 0.210*50 -> 10  ...  dava poradi mnoziny do ktere spadne, cislovano od 0
				z = (int)(rcon.p[j][3*i + 2] * gsi);

				(counts[x*gsi*gsi + y*gsi + z])++;
		}
	}
}

void point_density(std::vector<int> counts, voro::container_poly &rcon, int &gsi)
{
	// [out]	counts			vector of numbers of points in a given areas of lattice (of the same length as a vector of residuals)
	// [in]		con				the container with stored particles
	// [in]		gsi				the grid size (number of cubes in each direction)

	int i, j;
	int x, y, z;

	counts.resize(gsi*gsi*gsi);
	// loop over all generators -- mod its coordinates -- obtain j,k,l values -- inc appropriate value in the vector
	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box


			x = (int)(rcon.p[j][4 * i] * gsi);			// coordinates --> i,j,k values
			y = (int)(rcon.p[j][4 * i + 1] * gsi);		// (int) : 0.897*10 -> 8 , 0.210*50 -> 10  ...  dava poradi mnoziny do ktere spadne, cislovano od 0
			z = (int)(rcon.p[j][4 * i + 2] * gsi);

			(counts[x*gsi*gsi + y*gsi + z])++;
		}
	}
}



// fc ave_rad returns the average radius of particles in the container
double ave_rad(voro::container_poly &con)
{
	int i, j; 
	double arad = 0;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box
			arad = arad + con.p[j][4 * i + 3];
		}
	}

	arad = arad / con.total_particles();

	return arad;
}


void new_con(voro::container_poly &con, bool a)
{
	int nx, ny, nz;
	con.clear();

	pre_container_poly pcon(0, 1, 0, 1, 0, 1, true, true, true);  // true = periodic in given coordinate
																  
	if (a == false) { pcon.import("../data/Lag2.txt"); } else
	{ pcon.import("datacon.txt"); }
	
	pcon.guess_optimal(nx, ny, nz);  // guess								 
	pcon.setup(con);   // Set up the container class and import the particles from the pre-container 

}

int null_boxes(voro::container &con)
{
	int i, j;
	i = 0;
	for (j = 0; j < con.nxyz; j++) { // loop over boxes (con.nxyz = pocet boxu)
		
		if (con.co[j] == 0) { i++; } // con.co[j] = pocet castic v j-tem boxu

	}
	return i;
}


int null_boxes(voro::container_poly &con)
{
	int i, j;
	i = 0;
	for (j = 0; j < con.nxyz; j++) { // loop over boxes (con.nxyz = pocet boxu)

		if (con.co[j] == 0) { i++; } // con.co[j] = pocet castic v j-tem boxu

	}
	return i;
}


// fc empty_cells returns number of empty cells in the tessellation corresponding to the given container
int empty_cells(voro::container_poly &con)
{
	//	[in]	con		container with stored particles
	//	[out]	cells	the list (vector) of empty cell ids

	voronoicell c;
	int i, j, ni;
	ni = 0;
	//const long double PI = 3.141592653589793238L;
	//double constant = (4 / 3)*PI*pow(alfa, 3); // volume of the smallest cell with given hardcore parameter alfa
	bool cell;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == false) { ni++; } // the cell was not computed (it is empty)
			
		}
	}
	//std::cout << "(Non-computed (empty) cells: " << ni << ") ";

	return ni;
}

int nonempty_cells(voro::container_poly &con)
{
	//	[in]	con		container with stored particles
	//	[out]	cells	the list (vector) of empty cell ids

	voronoicell c;
	int i, j, ni;
	ni = 0;
	//const long double PI = 3.141592653589793238L;
	//double constant = (4 / 3)*PI*pow(alfa, 3); // volume of the smallest cell with given hardcore parameter alfa
	bool cell;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) { ni++; } // the cell was not computed (it is empty)

		}
	}
	//std::cout << "(Non-computed (empty) cells: " << ni << ") ";

	return ni;
}

void un_vertices(voro::container &rcon)
{
	unsigned int k;
	int i, j;
	int i0, i1, i2, i3, i4, i5, i6, i7, i8, i9;
	voronoicell c;  // bunka
	std::vector<int> ord;

	i0 = 0; i1 = 0; i2 = 0; i3 = 0; i4 = 0; i5 = 0; i6 = 0; i7 = 0; i8 = 0; i9 = 0;

	for (j = 0; j < rcon.nxyz; j++) { // loop over boxes
		for (i = 0; i < rcon.co[j]; i++) { // loop over particles in considered box

			rcon.compute_cell(c, j, i);
			c.vertex_orders(ord);

			for (k = 0; k < ord.size(); k++)
			{
				if (ord[k] == 0) { i0++; }
				if (ord[k] == 1) { i1++; }
				if (ord[k] == 2) { i2++; }
				if (ord[k] == 3) { i3++; }
				if (ord[k] == 4) { i4++; }
				if (ord[k] == 5) { i5++; }
				if (ord[k] == 6) { i6++; }
				if (ord[k] == 7) { i7++; }
				if (ord[k] == 8) { i8++; }
				if (ord[k] == 9) { i9++; }
			}

		}
	}
	std::cout << "Histogram of vertex orders (0-9): " << i0 << "  " << i1 << "  " << i2 << "  " << i3 << "  " << i4 << "  " << i5 << "  " << i6 << "  " << i7 << "  " << i8 << "  " << i9 << "\n";

}


// fce inter_boxes returns vector of box numbers of boxes intersecting the given sphere
bool inter_boxes(voro::container_poly &rcon, int ijk, int q, std::vector<int> boxno)
{
	//	[in]	con		examined container
	//	[in]	ijk,q	position of the particle - ijk is the number of the first box, q is the position within this box
	//	[out]	boxno	vector of box numbers of boxes intersecting the sphere

	boxno.clear();				//vycistit
	
	double rad, x_min, x_max, y_min, y_max, z_min, z_max;
	// box sizes: rcon.nx, rcon.ny, rcon.rz 

	rad = rcon.p[ijk][4 * q + 3];
	x_min = rcon.p[ijk][4 * q] - rad;
	x_max = rcon.p[ijk][4 * q] + rad;
	y_min = rcon.p[ijk][4 * q + 1] - rad;
	y_max = rcon.p[ijk][4 * q + 1] + rad;
	z_min = rcon.p[ijk][4 * q + 2] - rad;
	z_max = rcon.p[ijk][4 * q + 2] + rad;

	int x1, y1, z1, x2, y2, z2, x3, y3, z3;
	int step;

	x1 = (int)(x_min*rcon.nx);
	x2 = (int)(rcon.p[ijk][4 * q]*rcon.nx);
	x3 = (int)(x_max*rcon.nx)+1;
	y1 = (int)(y_min*rcon.ny);
	y2 = (int)(rcon.p[ijk][4 * q + 1] * rcon.nz);
	y3 = (int)(y_max*rcon.ny)+1;
	z1 = (int)(z_min*rcon.nz);
	z2 = (int)(rcon.p[ijk][4 * q + 2] * rcon.nz);
	z3 = (int)(z_max*rcon.nz)+1;

	int i, j, k;
	for (i = x1; i < x3; i++) {
		for (j = y1; j < y3; j++) {
			for (k = z1; k < z3; k++) {
				step = (i - 1)*rcon.ny*rcon.nz + (j - 1)*rcon.nz + k;
				boxno.push_back(step);
			}
		}
	}

	step = (x2 - 1)*rcon.ny*rcon.nz + (y2 - 1)*rcon.nz + z2;
	if (step == ijk) { return true; }
	else { std::cout << " Boxes do NOT fit! " << step << " vs " << ijk << " \n"; return false; }
}



void test_ed(voro::container &con)
{
	double x = uniform(0, 1);
	double y = uniform(0, 1);
	double z = uniform(0, 1);
	int id = con.total_particles() + 1;
	int ijk, q, j, k, nov;
	voronoicell_neighbor c;

	con.put(id, x, y, z);
	find_pos(ijk, q, id, &con);
	
	con.compute_cell(c, ijk, q);

	nov = 2 + c.number_of_edges() - c.number_of_faces(); // porovnat s rc1.p
	std::cout << c.p << "  vs  " << nov << "\n";

	for (j = 0; j < c.p; j++) {			// loop pres hrany: c.p = snad pocet vrcholu, c.nu - vektor radu vrcholu
		for (k = 0; k < c.nu[j]; k++) {
			std::cout << c.ed[j][k] << " ";
		}
		std::cout << "\n";
	}


}


void find05(voro::container_poly &con)
{
	int i, j, id;
	int cunt = 0;
	double x, y, z, r, xx, yy, zz;
	double h_min, h_max;
	voronoicell_neighbor c;
	std::vector<int> neigh;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			con.compute_cell(c, j, i);

			id = con.id[j][i];
			x = con.p[j][4 * i];
			y = con.p[j][4 * i + 1];
			z = con.p[j][4 * i + 2];
			r = con.p[j][4 * i + 3];

			h_fcs(c, x, y, z, x, y, z, h_max, h_min);

			if (h_max > 0.1) { 
				cunt++;  

				c.neighbors(neigh);
				c.centroid(xx, yy, zz);

				std::cout << h_max << ": " << id << " " << x << " " << y << " " << z << " " << r << " " << c.volume() << " " << c.number_of_faces() << " " << xx << " " << yy << " " << zz << " \n";
			}

		}
	}
	std::cout << " Total number of 05 cells: " << cunt << " \n";
}