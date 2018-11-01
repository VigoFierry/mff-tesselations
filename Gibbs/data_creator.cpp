#pragma warning (disable : 4996)

#include "Header.h"


// vytvoøení ètvercové sítì bodù v txt
void cube_net(double h)
{
	double i, j, k;
	const double alpha = h/2;
	int id = 1;

	FILE *f;
	f = fopen("../data/data.txt", "w");

	for (i = alpha; i <= 1; i += h) {
		for (j = alpha; j <= 1; j += h) {
			for (k = alpha; k <= 1; k += h) {

				fprintf(f, "%d %g %g %g \n", id,i,j,k);
				id += 1;
			}
		}
	}

	fclose(f);
}

// vytvoreni ctvercove mrize v txt, kde jsou navic doplneny vahy (defaultne vsechny 1)
void cube_rad_net(double h)
{
	double i, j, k;
	const double alpha = h / 2;
	int id = 1;
	double r = h/2;

	// C://Users/Vigo/Dokumenty/Visual Studio 2015/Projects/Gibbs-Voronoi/
	FILE *f;
	f = fopen("../data/data_L.txt", "w");

	for (i = alpha; i <= 1; i += h) {
		for (j = alpha; j <= 1; j += h) {
			for (k = alpha; k <= 1; k += h) {

				fprintf(f, "%d %g %g %g %g \n", id, i, j, k, r);
				id += 1;
			}
		}
	}

	fclose(f);
}

// kazdy bod mrize je nahodne posunut o malou vzdalenost
void random_net(double h)
{
	double i, j, k;
	double e1, e2, e3;
	const double alpha = h / 2;
	const double beta = h / 10;
	int id = 1;

	FILE *f;
	f = fopen("../data/datar.txt", "w");

	for (i = alpha; i <= 1; i += h) {
		for (j = alpha; j <= 1; j += h) {
			for (k = alpha; k <= 1; k += h) {

				e1 = uniform(-1, 1); e1 = e1*beta + i;
				e2 = uniform(-1, 1); e2 = e2*beta + j;
				e3 = uniform(-1, 1); e3 = e3*beta + k;

				fprintf(f, "%d %g %g %g \n", id, e1, e2, e3);
				id += 1;
			}
		}
	}

	fclose(f);
}

void random_rad_net(double h)
{
	double i, j, k;
	double e1, e2, e3, e4;
	const double alpha = h / 2;
	const double beta = h / 1000;
	int id = 1;
	double r = h/4;

	FILE *f;
	f = fopen("../data/datar_L.txt", "w");

	for (i = alpha; i <= 1; i += h) {
		for (j = alpha; j <= 1; j += h) {
			for (k = alpha; k <= 1; k += h) {

				e1 = uniform(-1, 1) * beta;
				e2 = uniform(-1, 1) * beta;
				e3 = uniform(-1, 1) * beta;
				e4 = uniform(-1, 1) * beta;

				fprintf(f, "%d %g %g %g %g \n", id, i + e1, j + e2, k + e3, r + e4);
				id += 1;
			}
		}
	}

	fclose(f);
}


// random container, only the total number of generators is specified apriori
void random_container(int n)
{
	double i;
	double x, y, z, r;
	int id = 1;

	FILE *f;
	f = fopen("../data/dataL_random.txt", "w");

	for (i = 0; i < n; i++) {
		x = uniform(0, 1);
		y = uniform(0, 1);
		z = uniform(0, 1);
		r = uniform(0, 1)*0.0625;

		fprintf(f, "%d %g %g %g %g \n", id, x, y, z, r);
		id += 1;
	}

	fclose(f);
}


void write_boxes(bool soubor, voro::container &con) {
	int i, j;

	if (soubor == 0) {
		for (j = 0; j < con.nxyz; j++) {
			std::cout << "BOX" << j << "  ";
			for (i = 0; i < con.co[j]; i++) {  // vypis castic v j-tem boxu
				std::cout << con.id[j][i] << " ";
			}
			std::cout << '\n';
		}
	}
	else {
		FILE *f;
		f = fopen("../data/boxes.txt", "w");
		for (j = 0; j < con.nxyz; j++) {
			fprintf(f, "BOX %d   ", j);
			for (i = 0; i < con.co[j]; i++) {
				fprintf(f, "%d ", con.id[j][i]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}
}

void write_boxes(bool soubor, voro::container_poly &con) {
	int i, j;

	if (soubor == 0) {
		for (j = 0; j < con.nxyz; j++) {
			std::cout << "BOX" << j << "  ";
			for (i = 0; i < con.co[j]; i++) {  // vypis castic v j-tem boxu
				std::cout << con.id[j][i] << " ";
			}
			std::cout << '\n';
		}
	}
	else {
		FILE *f;
		f = fopen("../data/boxes.txt", "w");
		for (j = 0; j < con.nxyz; j++) {
			fprintf(f, "BOX %d   ", j);
			for (i = 0; i < con.co[j]; i++) {
				fprintf(f, "%d ", con.id[j][i]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}

}

void write_container(voro::container &con) {
	int i, j;
	int citac = 0;

	FILE *f;
	f = fopen("data_sim.txt", "w");
	//f = fopen("../data/d200_1.txt", "w");
	if (f == NULL) { std::cout << "NELZE zapsat container \n"; }

	// vypise se id a souradnice castic podle boxu
	// 23.10.2017: ID castic z puvodnich dat nejsou aktualni - cislovani neodpovida poctu castic a postrada vyznam; nove proto budou generatory precislovany vzestupne od 1
	//				nejsou-li ID serazeny, nelze potom iterovat BDMA (resp vysledky jsou random)

	for (j = 0; j < con.nxyz; j++) {
		for (i = 0; i < con.co[j]; i++) {
			citac++;
			con.id[j][i] = citac;  // "online" precislovani
			fprintf(f, "%d %g %g %g \n", citac, con.p[j][3*i], con.p[j][3 * i + 1], con.p[j][3 * i + 2]);
		}
	}
	fclose(f);
}

void write_container(voro::container_poly &con) {
	int i, j;
	int citac = 0;

	FILE *f;
	f = fopen("dataconHIST.txt", "w");
	//f = fopen("../data/d200_1.txt", "w");
	if (f == NULL) { std::cout << "NELZE zapsat container \n"; }

	// vypise se id a souradnice castic podle boxu
	// 23.10.2017: ID castic z puvodnich dat nejsou aktualni - cislovani neodpovida poctu castic a postrada vyznam; nove proto budou generatory precislovany vzestupne od 1
	//				nejsou-li ID serazeny, nelze potom iterovat BDMA (resp vysledky jsou random)

	for (j = 0; j < con.nxyz; j++) {
		for (i = 0; i < con.co[j]; i++) {
			citac++;
			con.id[j][i] = citac;  // "online" precislovani
			fprintf(f, "%d %g %g %g %g \n", citac, con.p[j][4 * i], con.p[j][4 * i + 1], con.p[j][4 * i + 2], con.p[j][4 * i + 3]);
			// id, x, y, z, r
		}
	}
	fclose(f);
}

void transform() {
	int i, j;
	int nx, ny, nz;
	int citac = 0;

	voro::pre_container_poly pcon(0, 1, 0, 1, 0, 1, true, true, true);  // true = periodic in given coordinate
	pcon.import("C://Users/Vigo/Dokumenty/Visual Studio 2015/Projects/Gibbs-Voronoi/Gibbs-Voronoi/Simulations2018/data_sim.txt");
	pcon.guess_optimal(nx, ny, nz);  // guess
	voro::container_poly con(0, 1, 0, 1, 0, 1, nx, ny, nz, true, true, true, 8);
	pcon.setup(con);  // import  

	FILE *f;
	//f = fopen("data_simu.txt", "w");
	f = fopen("C://Users/Vigo/Dokumenty/Visual Studio 2015/Projects/Gibbs-Voronoi/Gibbs-Voronoi/Simulations2018/data_simu.txt", "w");
	//f = fopen("../data/d200_1.txt", "w");
	if (f == NULL) { std::cout << "NELZE zapsat container \n"; }

	// vypise se id a souradnice castic podle boxu
	// 23.10.2017: ID castic z puvodnich dat nejsou aktualni - cislovani neodpovida poctu castic a postrada vyznam; nove proto budou generatory precislovany vzestupne od 1
	//				nejsou-li ID serazeny, nelze potom iterovat BDMA (resp vysledky jsou random)

	for (j = 0; j < con.nxyz; j++) {
		for (i = 0; i < con.co[j]; i++) {
			citac++;
			con.id[j][i] = citac;  // "online" precislovani
			fprintf(f, "%d %g %g %g \n", citac, con.p[j][4 * i], con.p[j][4 * i + 1], con.p[j][4 * i + 2]);
			// id, x, y, z
		}
	}
	fclose(f);
}


/* Parameters for formating in fprintf:
d or i	Signed decimal integer	392
ld
u	Unsigned decimal integer	7235
o	Unsigned octal	610
x	Unsigned hexadecimal integer	7fa
X	Unsigned hexadecimal integer (uppercase)	7FA
f	Decimal floating point, lowercase	392.65
F	Decimal floating point, uppercase	392.65
e	Scientific notation (mantissa/exponent), lowercase	3.9265e+2
E	Scientific notation (mantissa/exponent), uppercase	3.9265E+2
g	Use the shortest representation: %e or %f	392.65
G	Use the shortest representation: %E or %F	392.65
a	Hexadecimal floating point, lowercase	-0xc.90fep-2
A	Hexadecimal floating point, uppercase	-0XC.90FEP-2
c	Character	a
s	String of characters	sample
p	Pointer address	b8000000
*/

/*
newline				\n
horizontal tab		\t
vertical tab		\v
backspace			\b
carriage return		\r
form feed			\f
alert				\a
backslash			\\
question mark		? or \?
single quote		\'
double quote		\"
the null character	\0
...
*/