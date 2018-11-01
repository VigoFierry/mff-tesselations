#include<iostream>
#include<cmath>
#include<iomanip>

#include "Header.h"


using namespace std;


// firstly define the function to be solved; 
double f(double x, std::vector<int> &ra, std::vector<double> &rb, std::vector<double> &rc, int &rk, int &ry)    //define the function here, ie give the equation
{
	// [in]		x		function variable
	// [in]		a,b,c,k	parameters of the function
	// [in]		y		y=0 ... MPLE ; y=1 ... "quick" estimate

	int i;
	unsigned int j;
	double fv = 0;  
	double k = 0;
	double theta0 = 1, theta1 = 1;			//zname konstanty

	if (ry == 0) {  // fc type 1 (MPLE)
		// fc is the derivative of PLL function w.r.t. theta with expression of zet obtained from the derivation of PLL w.r.t. zet
		for (i = 0; i < rk; i++) { k = k + rb[i]; }
		k = (k / rk);
		for (j = 0; j < rc.size()/3; j++) {

			if (ra[j] == 1) {
			 // fv = fv + ((exp(-x*rc[j])) * (rc[j] - k)); // ORIGINAL for one potential
				fv = fv + ((exp(theta0*rc[3*j]+theta1*rc[3*j+1]+x*rc[3*j+2])) * (-rc[3*j+2] - k)); // !!! POZOR na znaminka
			}
		}

		return fv;
	}
	else {	// fc type 2 (QUICK)  ... ORIGINAL OK, uprava ???
		for (i = 0; i < rk; i++) { k = k + exp(x*rb[i]); }
		for (j = 0; j < rc.size(); j++) {

			if (ra[j] == 1) {
			//	fv = fv + exp(-x*rc[j]); // ORIGINAL for one potential
				fv = fv + exp(theta0*rc[3 * j] + theta1*rc[3 * j + 1] + x*rc[3 * j + 2]);
			}
		}
		fv = fv*k - ry*rk;

		return fv;
	}
}

// for the Newton-Raphson method we will need its first derivative too
double fprime(double x, std::vector<int> &ra, std::vector<double> &rb, std::vector<double> &rc, int &rk, int &ry)
{
	int i;
	unsigned int j;
	double fvv = 0;        
	double k = 0;
	double l = 0;
	double m = 0;

	double theta0 = 1, theta1 = 1;			//zname konstanty
	
	if (ry == 0) {  // fc type 1 (MPLE)
		for (i = 0; i < rk; i++) { k = k + rb[i]; }
		k = k / rk;
		for (j = 0; j < rc.size(); j++) {

			if (ra[j] == 1) {
				//fvv = fvv + ((exp(-(x*rc[j]))) * ((rc[j])*(k - rc[j]))); // ORIGINAL for one potential
				fvv = fvv + (exp(theta0*rc[3 * j] + theta1*rc[3 * j + 1] + x*rc[3 * j + 2])) * ((-rc[3*j+2])*(k + rc[3*j+2]));
			}
		}

		return fvv;
	}
	else {	// fc type 2 (QUICK) ... ORIGINAL OK, uprava ???
		for (i = 0; i < rk; i++) { 
			k = k + ((exp(x*rb[i]))*rb[i]);
			m = m + exp(x*rb[i]);
		}
		for (j = 0; j < rc.size(); j++) {

			if (ra[j] == 1) {
				fvv = fvv - ((exp(-(x*rc[j])))*rc[j]);
				l = l + exp(-(x*rc[j]));
			}
		}
		fvv = ((fvv*m + l*k)/ry);

		return fvv;
	}
}

// methods for finding the numeric solution of given algebraic/transcendental equation

// 1) bisection method
double bisection(std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int &k, int y)
{
cout.precision(4);        //set the precision
cout.setf(ios::fixed);
double a, b, c, e, fa, fb, fc;    //declare some needed variables
a:cout << "Enter the initial guesses:\na=";    //Enter the value of a(set a label('a:') for later use with goto)
cin >> a;
cout << "\nb=";            //Enter the value of b
cin >> b;
// cout << "\nEnter the degree of accuracy desired" << endl;    //Enter the accuracy
// cin >> e;                //e stands for  accuracy
e = 0.01;
if (f(a, va, vb, vc, k, y)*f(b, va, vb, vc, k, y)>0)        //Check if a root exists between a and b
{                //If f(a)*f(b)>0 then the root does not exist between a and b
	cout << "Please enter a different intial guess" << endl;
	goto a;            //go back to 'a' ie 17 and ask for different values of a and b
}
else                //else a root exists between a and b
{
	while (fabs(a - b) >= e)        /*if the mod of a-b is greater than the accuracy desired keep bisecting the interval*/
	{
		c = (a + b) / 2.0;        //bisect the interval and find the value of c
		fa = f(a, va, vb, vc, k, y);
		fb = f(b, va, vb, vc, k, y);
		fc = f(c, va, vb, vc, k, y);
		cout << "a=" << a << "     " << "b=" << b << "     " << "c=" << c << "      fc=" << fc << endl;/*print the values of a,b,c and fc  after each iteration*/
		if (fc == 0)        //if f(c)=0, that means we have found the root of the equation
		{
			cout << "The root of the equation is " << c;    /*print the root of the equation and break out of the loop*/
			break;
		}

		if (fa*fc>0)    //if f(a)xf(c)>0, that means no root exist between a and c 
		{
			a = c;    /*hence make a=c, ie make c the starting point of the interval and b the end point*/
		}
		else if (fa*fc<0)
		{
			b = c;    /*this means that a root exist between a and c therfore make c the end point of the interval*/
		}


	}
}            //The loop ends when the difference between a and b becomes less than the desired accuracy ie now the value stored in 'c' can be called the approximate root of the equation         
cout << "The root of the equation is " << c;    //print the root    
return c;
}

// 2) Secant Method for finding the roots of an equation
double secant(std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int &k, int y)
{
	cout.precision(4);
	cout.setf(ios::fixed);        //set the precision of the output
	double a, b, c, e;
	cout << "Enter the initial guess\na=";
	cin >> b;
	cout << "b=\n";                //take an intial guess
	cin >> c;
	// cout << "Enter the degree of accuracy\n";
	// cin >> e;                    //take the desired accuracy
	e = 0.01;
	do
	{
		a = b;
		b = c;                //make b equal to the last calculated value of c
		c = b - (b - a) / (f(b, va, vb, vc, k, y) - f(a, va, vb, vc, k, y))*f(b, va, vb, vc, k, y);    //calculate c
		if (f(c, va, vb, vc, k, y) == 0)
		{
			cout << "\nThe root of the equation is " << c;    //print the root
			return 0;
		}
	} while (abs(c - b) >= e);            //check if the error is greater than the desired accuracy
	cout << "\nThe root of the equation is " << c;    //print the root
	return c;
}

// 3) Newton-Raphson Method
double NR(std::vector<int> &va, std::vector<double> &vb, std::vector<double> &vc, int &k, int y)
{

	double x, x1, e, fx, fx1;
	cout.precision(4);        //set the precision
	cout.setf(ios::fixed);
	cout << "NR: Enter the initial guess\n";    //take an intial guess
	// cin >> x1; // set the initial guess manually
	x1 = 0; // set the initial guess defaultly
	// cout << "Enter desired accuracy\n";    //take the desired accuracy
	// cin >> e;
	e = 0.005;
	//	fx = f(x,va,vc,k);
	//	fx1 = fprime(x,va,vc);
	cout << "x{i}" << "    " << "x{i+1}" << "        " << "|x{i+1}-x{i}|" << endl;

	do
	{
		x = x1;                /*make x equal to the last calculated value of x1*/
		fx = f(x,va,vb,vc,k,y);            //simplifying f(x)to fx
		fx1 = fprime(x,va,vb,vc,k,y);            //simplifying fprime(x) to fx1
		x1 = x - (fx / fx1);            /*calculate x{1} from x, fx and fx1*/
		//cout << x << "     " << x1 << "           " << abs(x1 - x) << endl;
	} while (fabs(x1 - x) >= e);            /*if |x{i+1}-x{i}| remains greater than the desired accuracy, continue the loop*/
	cout << "The root of the equation is " << x1 << endl;
	return x1;
}

