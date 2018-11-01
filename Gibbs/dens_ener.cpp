#include "Header.h"


using namespace voro;


// class HISTOGRAM

// histogram = vector; [0] ... beginning value, [1] ... step size, the rest contains number of occurencies
//				it should have a standard length; 
//				the range should be that one observed from data (values outside the range are not important, everything we need to know is its number)
//   ... class histogram



// constructor:
histogram::histogram(double sval, double tep, std::vector<int> v) {

	sp = sval;
	step = tep;
	oc = v;

	int max = 0;
	int sum = 0;
	for (int i = 0; i < v.size(); i++) {
		if (v[i] > max) { max = v[i]; }
		sum = sum + v[i];
	}

	ocm = max;
	so = sum;
	noc = 0;
}

void histogram::read_hist_vol() 
{
	
	std::ifstream infile;
	//std::ifstream infile("C:/Random tesselations/Data/data Laguerre/Krill data/hist_rad.txt");
	//infile.open("C:/Random tesselations/Data/data Laguerre/Krill data/hist_rad.txt");
	//infile.open("C:/Random tesselations/Data/data Laguerre/Krill data/hist_vol.txt");  // 113 sloupcu
	// !!! pri nacitani jineho histogramu je potreba zmenit pocet sloupcu v create_hist_double
	infile.open("../data/hist_vol.txt");

	if (!infile) {
		std::cout << "ERROR: (read_hist) NELZE cist \n";
		//exit(1);   // call system to stop
	}

	int nom, val;
	infile >> sp;
	infile >> step;
	infile >> ocm;
	infile >> nom;

	//std::cout << sp << " - hahaha \n";
	//std::cout << step << " - hahaha \n";
	//std::cout << ocm << " - hahaha \n";
	//std::cout << nom << " - hahaha \n";

	int sum = 0;

	oc.clear();
	for (int i = 0; i < nom; i++) {		// read vector (ints)
		infile >> val;
		oc.push_back(val);
		sum = sum + val;
	}

	so = sum;
	noc = 0;

	infile.close();
}


void histogram::read_hist_nof()
{

	std::ifstream infile;
	//std::ifstream infile("C:/Random tesselations/Data/data Laguerre/Krill data/hist_rad.txt");
	//infile.open("C:/Random tesselations/Data/data Laguerre/Krill data/hist_rad.txt");
	//infile.open("C:/Random tesselations/Data/data Laguerre/Krill data/hist_nof.txt");	// 32 sloupcu
	// !!! pri nacitani jineho histogramu je potreba zmenit pocet sloupcu v create_hist_int
	infile.open("../data/hist_nof.txt");

	if (!infile) {
		std::cout << "ERROR: (read_hist) NELZE cist \n";
		//exit(1);   // call system to stop
	}

	int nom, val;
	infile >> sp;
	infile >> step;
	infile >> ocm;
	infile >> nom;

	//std::cout << sp << " - hahaha \n";
	//std::cout << step << " - hahaha \n";
	//std::cout << ocm << " - hahaha \n";
	//std::cout << nom << " - hahaha \n";

	int sum = 0;

	oc.clear();
	for (int i = 0; i < nom; i++) {		// read vector (ints)
		infile >> val;
		oc.push_back(val);
		sum = sum + val;
	}

	so = sum;
	noc = 0;

	infile.close();
}

int histogram::hist_value(double rc)
{
	if (rc < sp) { return 0; }

	int size = oc.size();
	if ((sp + size*step) <= rc) { return 0; }

	for (int i = 0; i < size; i++) {
		if (sp + i*step <= rc && rc < sp + (i + 1)*step) { return oc[i]; }
	}

	return 0; // predchozi pripady jsou METE
}

void histogram::create_hist_int(voro::container_poly &con)
{
	// [in]		con		container

	// jako prvni je potreba specifikovat geometrickou charakteristiku, napr. volume, number of neighbours, ...

	// dale je potreba specifikovat parametry histogramu (step value, beginning value)
	// pozn.: beginning value by mela byt defaultni (asi 0), prtz kdyz budu chtit pouzit tuto fci na pocatecni konfiguraci, musim pocitat
	//		s tim, ze pocatecni konfigurace muze byt libovolna

	// ??? spocitat nejdriv statistiky vsech bunek a pak az prevest na histogram (vyhodnejsi asi pro objem - neznam rozpeti hodnot),
	//	nebo aktualizovat cetnost v histogramu ihned po spocteni statistiky kazde bunky (tj predefinovat histogram predem, vyhodou pro nof) 

	int i, j;
	voronoicell c;
	bool cell;
	std::vector<double> vols;
	int nof, size;
	int max = 0;

	// values given by the histogram from real data:
	sp = 4;	// minimal value for number of faces / vertices in 3D
	step = 1; // step for integers
	size = 32;
	noc = 0;

	oc.clear();
	oc.resize(size);

	int n = 0;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) { 
				// compute statistic
				nof = c.number_of_faces();
				//if (nof > max) { max = nof; }
				if (nof > sp + size - 1) { noc++; } // increase histogram	
				else { oc[nof - 4]++; }
				//if (nof == 5) { n++; }
				n++;
			} 

		}
	}
	
	//oc.resize(max - 4);
	so = n;

	max = 0;
	for ( i = 0; i < oc.size(); i++) {
		if (oc[i] > max) { max = oc[i]; }
	}
	
	ocm = max;
}


void histogram::create_hist_double(voro::container_poly &con)
{
	// [in]		con		container

	// jako prvni je potreba specifikovat geometrickou charakteristiku, napr. volume, number of neighbours, ...

	// dale je potreba specifikovat parametry histogramu (step value, beginning value)
	// pozn.: beginning value by mela byt defaultni (asi 0), prtz kdyz budu chtit pouzit tuto fci na pocatecni konfiguraci, musim pocitat
	//		s tim, ze pocatecni konfigurace muze byt libovolna

	// ??? spocitat nejdriv statistiky vsech bunek a pak az prevest na histogram (vyhodnejsi asi pro objem - neznam rozpeti hodnot),
	//	nebo aktualizovat cetnost v histogramu ihned po spocteni statistiky kazde bunky (tj predefinovat histogram predem, vyhodou pro nof) 

	int i, j;
	voronoicell c;
	bool cell;
	std::vector<double> vols;
	double vol;
	int nof, size;
	//int cou = 50; // number of histogram boxes  - must be computed from step
	double max = 0, min = 1000000000;

	vols.clear();
	vols.reserve(con.total_particles());

	// !!!!!!!!!!!!! step is always GIVEN by histogram coming from real data
	step = 0.00005; 
	// beginning value, step and oc.size have to be the same as in the histogram coming from real data (otherwise we loose the ability to compare)
	sp = 0; 
	size = 113;

	oc.resize(size);

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {
				// compute statistic
				vol = c.volume();
				if (vol > max) { max = vol; }
				if (vol < min) { min = vol; }
				vols.push_back(vol);
			}

		}
	}

	double cou = floor(min / step);

// -	sp = cou*step; // beginning value is the biggest multiple of step smaller than min
	// !!!!!!!!!!!!!!!! step = ((max - min) + (max - min) / (4 * cou)) / cou;

	// beginning value, step and oc.size have to be the same as in the histogram coming from the real data
	
	// number of boxes in histogram 
// -	cou = ceil(max / step) - cou;
// -	oc.resize((int)cou);

	so = vols.size();
	noc = 0;

	for (i = 0; i < so; i++) {
		//j = (int)floor((vols[i] - sp) / step);
		j = static_cast<int>(floor((vols[i] - sp) / step));
		if (j > -1 && j < size) { oc[j]++; }
		else { noc++; }
	}

	max = 0;
	for (i = 0; i < oc.size(); i++) {
		if (oc[i] > max) { max = oc[i]; }
	}

	ocm = max;
}


// funkce zvysujici/snizujici cetnost pro danou hodnotu
void histogram::hist_act(double val, bool op) { 

	int size = oc.size();
	
	if (val < sp) { 
		if (op == 1) { noc++; so++; } else { noc--; so--; }
		return; }
	if (val > sp + size*step) { 
		if (op == 1) { noc++; so++; } else { noc--; so--; }
		return; }

	for (int i = 0; i < size; i++) {
		if (sp + i*step <= val && val < sp + (i + 1)*step) { 
			//std::cout << val << " " << i << " \n";
			if (op == 1) { oc[i]++; so++;  if (oc[i] == ocm) { ocm++; } }
			else { 
				oc[i]--; 
				so--;
				if (oc[i] == ocm) { 
					ocm--;		// muze existovat j neq i tak, ze oc[j] = ocm, pak by ocm zustalo stejne
					for (int j = 0; j < oc.size(); j++) {
						if (oc[j] > ocm) { ocm = oc[j];  break; } // muze se zvysit jen jednou o 1
					}
				}  
			}
		}
	}

	// zmenou cetnosti se muze zmenit i hodnota ocm 
	// ocm je jedna z moznosti jak histogram normovat (dalsi moznosti je napr soucet vsech cetnosti)

	// --> ocm se aktualizuje hur nez prosty soucet cetnosti !!!!!

}

/*
// fc returns disrepancy between two histograms (absolute one)
double hist_dis(histogram &hist1, histogram &hist2) {

	// [in] hist1
	// [in] hist2

	// assumption: the step of both histograms is equal and both beginning values are multiples of this step
	// (step is always given by the histogram from real data)
	if (hist1.step == hist2.step) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different step. \n";  return 0; }

	bool frst, lst;
	int m, l;
	int dist = 0, i;
	int s1 = hist1.oc.size(), s2 = hist2.oc.size();
	if (hist1.sp < hist2.sp) { 
		frst = 0;
		m = (int)((hist2.sp - hist1.sp) / hist1.step);
		if (s1 < (s2 + m)) {
			lst = 1;
			l = s1 - m;
			for (i = 0; i < s2 - l; i++) { dist = dist + hist2.oc[i+l]; }
		} else {
			lst = 0;
			l = s2;
			for (i = 0; i < s1 - l - m ; i++) { dist = dist + hist1.oc[i+l+m]; }
		}

		for (i = 0; i < m; i++) { dist = dist + hist1.oc[i]; }
		for (i = 0; i < l; i++) { dist = dist + abs(hist1.oc[i+m]-hist2.oc[i]); }
		
	} else { 
		frst = 1; 
		m = (hist1.sp - hist2.sp) / hist1.step;
		if (s2 < (s1 + m)) {
			lst = 0;
			l = s2 - m;
			for (i = 0; i < s1 - l; i++) { dist = dist + hist1.oc[i + l]; }
		}
		else {
			lst = 1;
			l = s1;
			for (i = 0; i < s2 - l - m; i++) { dist = dist + hist2.oc[i + m + l]; }
		}

		for (i = 0; i < m; i++) { dist = dist + hist2.oc[i]; }
		for (i = 0; i < l; i++) { dist = dist + abs(hist2.oc[i + m] - hist1.oc[i]); }
	}

	return (double)dist;
}


// fc returns disrepancy between two histograms (proportional)
double hist_disp(histogram &hist1, histogram &hist2) {

	// [in] hist1
	// [in] hist2

	// assumption: the step of both histograms is equal and both beginning values are multiples of this step
	// (step is always given by the histogram from real data)
	if (hist1.step == hist2.step) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different step. \n";  return 0; }

	bool frst, lst;
	int m, l, i;
	double dist = 0;
	int s1 = hist1.oc.size(), s2 = hist2.oc.size();
	if (hist1.sp < hist2.sp) {
		frst = 0;
		m = (int)((hist2.sp - hist1.sp) / hist1.step);
		std::cout << m << "\n";
		if (s1 < (s2 + m)) {
			lst = 1;
			l = s1 - m;
			for (i = 0; i < s2 - l; i++) { dist = dist + (double)hist2.oc[i + l] / hist2.so; }
		}
		else {
			lst = 0;
			l = s2;
			for (i = 0; i < s1 - l - m; i++) { dist = dist + (double)hist1.oc[i + m + l + 1] / hist1.so; }
		}

		for (i = 0; i < m; i++) { dist = dist + (double)hist1.oc[i] / hist1.so; }
		for (i = 0; i < l; i++) { dist = dist + abs((double)hist1.oc[i + m]/hist1.so - (double)hist2.oc[i]/hist2.so); }

	}
	else {
		frst = 1;
		m = (hist1.sp - hist2.sp) / hist1.step;
		std::cout << m << "\n";
		if (s2 < (s1 + m)) {
			lst = 0;
			l = s2 - m;
			for (i = 0; i < s1 - l; i++) { dist = dist + (double)hist1.oc[i + l]/hist1.so; }
		}
		else {
			lst = 1;
			l = s1;
			for (i = 0; i < s2 - l - m; i++) { dist = dist + (double)hist2.oc[i + m + l]/hist2.so; }
		}

		for (i = 0; i < m; i++) { dist = dist + (double)hist2.oc[i]/hist2.so; }
		for (i = 0; i < l; i++) { dist = dist + abs((double)hist2.oc[i + m]/hist2.so - (double)hist1.oc[i]/hist1.so); }
	}

	return dist;
}*/

// fc returns disrepancy between two histograms (absolute)
double hist_dis(histogram &hist1, histogram &hist2) {

	// [in] hist1	comparable histogram
	// [in] hist2	original histogram from data

	// assumption: the step, beginning value and oc.size of both histograms is equal 
	// (step, beginning value and oc.size are always given by the histogram from real data)
	if (hist1.step == hist2.step) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different step. \n";  return 0; }

	if (hist1.sp == hist2.sp) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different beginning value. \n";  return 0; }

	if (hist1.oc.size() == hist2.oc.size()) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different oc.size. \n";  return 0; }

	if (hist2.noc == 0) {}
	else { std::cout << "ERROR: (hist_dis) histogram is not original. \n";  return 0; }

	int dist = 0;
	int l = hist1.oc.size();

//	for (int ii = 0; ii < l; ii++) { std::cout << hist1.oc[ii] << " "; } 	std::cout << "\n";
//	for (int ii = 0; ii < l; ii++) { std::cout << hist2.oc[ii] << " "; } 	std::cout << "\n";

	for (int i = 0; i < l; i++) { dist = dist + abs(hist1.oc[i] - hist2.oc[i]); }

	//dist = dist + (int)hist1.noc;
	dist = dist + static_cast<int>(hist1.noc);

	//return (double)dist;
	return static_cast<double>(dist);

}

// fc returns disrepancy between two histograms (proportional)
double hist_disp(histogram &hist1, histogram &hist2) {

	// [in] hist1	comparable histogram
	// [in] hist2	original histogram from data

	// assumption: the step, beginning value and oc.size of both histograms is equal 
	// (step, beginning value and oc.size are always given by the histogram from real data)
	if (hist1.step == hist2.step) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different step. \n";  return 0; }

	if (hist1.sp == hist2.sp) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different beginning value. \n";  return 0; }

	int l = hist1.oc.size();
	int l2 = hist2.oc.size();

	if (l == l2) {}
	else { std::cout << "ERROR: (hist_dis) histograms with different oc.size. \n";  return 0; }

	if (hist2.noc == 0) {}
	else { std::cout << "ERROR: (hist_dis) histogram is not original. \n";  return 0; }

	double dist = 0;
	

	//for (int i = 0; i < l; i++) { dist = dist + abs((double)hist1.oc[i] / hist1.so - (double)hist2.oc[i] / hist2.so); }
	for (int i = 0; i < l; i++) { dist = dist + abs_val(static_cast<double>(hist1.oc[i]) / hist1.so - static_cast<double>(hist2.oc[i]) / hist2.so); }


	dist = dist + hist1.noc / hist1.so; 

	return dist;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// class CON_INFO

// constructor:
con_info::con_info(int n, std::vector<double> v, std::vector<double> w, histogram hist) {

	tp_bef = n;
	tp_aft = n;
	mean_bef = v;
	mean_aft = v;
	var_bef = w;
	var_aft = w;
	hist_bef = hist;
	hist_aft = hist;
	empty.clear();
}

void con_info::get_meansum(voro::container_poly &con) {

	bool cell;
	voronoicell c;
	int t;
	std::vector<double> xn;

	xn.push_back(0);
	xn.push_back(0);

	for (int j = 0; j < con.nxyz; j++) { // loop over boxes
		for (int i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {
				xn[0] = xn[0] + c.number_of_faces();		// nof
				xn[1] = xn[1] + c.volume();					// volume
				// ...
			} // the cell was computed 

		}
	}

	mean_bef.clear();
	mean_aft.clear();

	//tp = nonempty_cells(con);

	for (int i = 0; i < xn.size(); i++) {
		mean_bef.push_back(xn[i]);
		mean_aft.push_back(xn[i]);
	}
}


void con_info::get_varsum(voro::container_poly &con) {

	bool cell;
	voronoicell c;
	std::vector<double> xn;

	xn.push_back(0);

	for (int j = 0; j < con.nxyz; j++) { // loop over boxes
		for (int i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {
				xn[0] = xn[0] + pow(c.number_of_faces() - mean_bef[0]/tp_bef,2);		// nof
				// xn[1] = xn[1] + pow(c.volume() - mean_bef[1]/tp_bef,2);				// volume
															// ...
			} // the cell was computed 

		}
	}

	var_bef.clear();
	var_aft.clear();

	//tp = nonempty_cells(con);

	for (int i = 0; i < xn.size(); i++) {
		var_bef.push_back(xn[i]);
		var_aft.push_back(xn[i]);
	}
}


void con_info::varsum(voro::container_poly &con) {

	bool cell;
	voronoicell c;
	std::vector<double> xn;

	xn.push_back(0);

	for (int j = 0; j < con.nxyz; j++) { // loop over boxes
		for (int i = 0; i < con.co[j]; i++) { // loop over particles in considered box

			cell = con.compute_cell(c, j, i);
			if (cell == true) {
				xn[0] = xn[0] + pow(c.number_of_faces() - mean_bef[0] / tp_bef, 2);		// nof
																						// xn[1] = xn[1] + pow(c.volume() - mean_bef[1]/tp_bef,2);				// volume
																						// ...
			} // the cell was computed 

		}
	}

	var_aft.clear();

	//tp = nonempty_cells(con);

	for (int i = 0; i < xn.size(); i++) {
		var_aft.push_back(xn[i]);
	}
}


void delete_empty(voro::container_poly &con) {

	bool cell;
	voronoicell c;
	int ci, k, i, j;
	//int id, ijk, q;
	int citac;


	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		ci = 0; k = con.co[j]; citac = 0;
		//std::cout << k << " con.co[j] \n";  /////////////////////////////////////////////////////////////////////////////////
		for (i = 0; i < k; i++) { // loop over particles in considered box

			ci = i - citac;

			//id = con.id[j][ci];
			//std::cout << j << " " << ci << " " << i << " \n"; ///////////////////////////////////////////////////////////////////////////
			//std::cout << id << " ";  /////////////////////////////////////////////////////////////////////////////////
			//find_pos(ijk, q, id, &con);	// find position of the particle
			//std::cout << ijk << " " << q << " \n"; ///////////////////////////////////////////////////////////////////////////

			cell = con.compute_cell(c, j, ci);
			if (cell == false) {
				// delete it for free
				erase(j, ci, &con);
				citac = citac + 1;
			//	std::cout <<  citac << " deleted \n";  /////////////////////////////////////////////////////////////////////////////////

				
			} // the cell is empty 
			
		}
	}

}


