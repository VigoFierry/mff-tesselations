// Generating random variables
// Date     : March 19th 2017

/* #include <iostream>
#include <random>  // package for generating random variables - http://en.cppreference.com/w/cpp/numeric/random or http://itscompiling.eu/2016/04/11/generating-random-numbers-cpp/
*/

#include "Header.h"
// prototypes of functions:
/*
double rnd();
double uniform(int, int);
int uniform_int(int a, int b);
double normal(int, int);
double poisson(int);
*/


/* ukazka pouziti:
int main()
{
int x;

std::cout << uniform(0, 1);
std::cout << '\n';

std::cout << normal(0, 1);
std::cout << '\n';

std::cout << poisson(1);
std::cout << '\n';

std::cin >> x;


return 0;
}	
	*/


	
double rnd()  // This function returns a random floating point number between 0 and 1 - NEPOUZIVAT RADSI
{ 
	srand(time(NULL));
	return double(rand()) / RAND_MAX; 
} 

double uniform(int a, int b)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	std::uniform_real_distribution<> dis(a, b);
	return dis(gen);

}

int uniform_int(int a, int b)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	std::uniform_int_distribution<> dis(a, b);
	return dis(gen);

}

double normal(double mu, double sigma) // jake typy je mozne pouzit pro parametry???????
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	std::normal_distribution<> d(mu, sigma);
	return d(gen);
}

double poisson(int lambda)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	std::poisson_distribution<> d(lambda);
	return d(gen);

}

double gamma(double alfa, double beta)
{
	double rate = 1 / beta;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	std::gamma_distribution<> d(alfa,rate);
	return d(gen);

}

double triangle(double a, double b, double c)		// triangle distribution
{
	double F = (c - a) / (b - a);
	double r = uniform(0, 1);						// using inverse tranformation method
	if ((a < c) & (c < b)) { 
		if (r < F) { return (a + sqrt(r*(b - a)*(c - a))); }
		else { return (b - sqrt((1 - r)*(b - a)*(b - c))); }
	}
	else { std::cout << "Please enter correct parameters (i.e. a < c < b) \n"; return -1; }
}

double exponential(double lambda)
{
	double r = uniform(0, 1);						// using inverse tranformation method

	return (- log(r)/lambda);
}

// NOTE: parametry INT mnohde stacit nebudou !!!