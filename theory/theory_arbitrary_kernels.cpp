/******************************************
The code for theoretical predictions in the exact combinatorial approach to coagulation - arbitrary kernels
Author: Michal Lepek
Date: 12 Feb 2024
If using this code for your research, please refer to the main work:
"A review on combinatorial approach to aggregation"
by M. Lepek, A. Fronczak & P. Fronczak, submitted to Physica D: Nonlinear Phenomena.
******************************************/

/// *********** Useful note:
// To calculate theoretical predictions for a given kernel, you must:
// 1. Set appropriate Kernel in the expression for {xn} - define or choose from the predefined equations that I used before (lines 93-115).
// 2. All further expressions, i.e. for the average cluster number of a given size, <n_g>, and its standard deviation, are ready to use and no need to modify them.
// 3. Remember to set proper output filenames (lines 187-188).
// 4. And of course set the system size M and the time of observation t (lines 34-35).

#include <iostream>
#include <gmpxx.h> // mpfr library
#include <vector>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <mpreal.h> // mpfrc++ library

using namespace std;
using mpfr::mpreal;

typedef vector<int> IVector;

mpreal Factorial(size_t n);
mpreal BinomialCoefficient(int n, int k);
mpreal Bell_g_k_specific_xn(int g, int k);

const int M = 250; // number of all monomers in the system
const int t = 224;  // time when the system is observed

const int k = M - t; // number of clusters observed in the time t (this code works for t < M !!)
const int max_g = M - k + 1; // size of the biggest cluster that can be potentially observed

// The 4 objects below are used for storing partial results during calculations of Bell polynomials.
// It is for optimization purposes: the function Bell_g_k_specific_xn() can call the same Bell many times.
// If this specific Bell polynomial was calculated before, it is taken from the memory, not calculated again.
// Without this feature, computing time would be unacceptable long.

// table with information if the Bell polynomial has been calculated before
int done_Bells[M+1][M+1];
// table with the Bell polynomial values
mpreal calculated_Bells[M+1][M+1];

// table with information if given Bell coefficient has been calculated before
int done_coeffs[M+1+1+1][M+1];
// table with the Bell coefficient values
mpreal calculated_coeffs[M+1+1+1][M+1];

// This vector will hold values of the mathematical sequence {xn} which is the argument to Bell polynomial
vector<mpreal> x;

//************************* main() ******************************
int main()
{
    	// heuristic value of precision which is enough for these calculations
	mpfr::mpreal::set_default_prec(50000);

	// Here, we calculate the mathematical sequence {xn} which is the argument to Bell polynomial.
	// For a given kernel, you must define {xn}.

	cout << "Calculating {xn}" << endl;

	// This vector holds values of the
	x = vector<mpreal>(0);

	// The first element of the x vector does not matter,
	// it is added only for starting useful values from x.at(1) instead x.at(0)
	x.push_back(1);
	
	// This element is for g=1
	x.push_back(1);
	
	// ... so here I start from i=2
	for(long int i=2; i<=M; ++i)
	{
		cout << "x" << i << endl;

		// Below, I calculate the sequence {x_n} by using a general recursive expression
		
		// Masses of the two merging sub-clusters are k and (i-k) !
		
		mpreal sum = 0;
		
		for(int k=1; k<=i-1; k++)
		{
		
			// Expression for the constant kernel
			sum += Factorial(i) / Factorial(k)  *  Factorial(i-2) / Factorial(i-k)  *  x.at(k) / Factorial(k-1)  *   x.at(i-k) / Factorial(i-k-1) * 1;
			
			// Expression for the product kernel
			//sum += Factorial(i) / Factorial(k)  *  Factorial(i-2) / Factorial(i-k)  *  x.at(k) / Factorial(k-1)  *   x.at(i-k) / Factorial(i-k-1)  * k*(i-k);
			
			// Expression for the additive kernel
			//sum += Factorial(i) / Factorial(k)  *  Factorial(i-2) / Factorial(i-k)  *  x.at(k) / Factorial(k-1)  *   x.at(i-k) / Factorial(i-k-1)  * i;
			
			// Expression for the linear-chain kernel
			//sum += Factorial(i) / Factorial(k)  *  Factorial(i-2) / Factorial(i-k)  *  x.at(k) / Factorial(k-1)  *   x.at(i-k) / Factorial(i-k-1)  * (1./k + 1./(i-k));
			
			// Brownian motion - continuum regime
			//sum += Factorial(i) / Factorial(k)  *  Factorial(i-2) / Factorial(i-k)  *  x.at(k) / Factorial(k-1)  *   x.at(i-k) / Factorial(i-k-1)  * ( mpfr::cbrt(k) + mpfr::cbrt(i-k) ) * ( mpfr::cbrt(1./k) + mpfr::cbrt(1./(i-k)) ); 
			
			// Brownian motion - free molecular regime
			//sum += Factorial(i) / Factorial(k)  *  Factorial(i-2) / Factorial(i-k)  *  x.at(k) / Factorial(k-1)  *   x.at(i-k) / Factorial(i-k-1)  * (pow(k, 1./3.)+pow((i-k),1./3.))*(pow(k, 1./3.)+pow((i-k),1./3.)) * sqrt(1./k + 1./(i-k));
		
			// Planetesimal kernel for creation of planetary embryos
			//double mju = 60;
			//mpreal i_mpreal = i;
			//mpreal K = ( cbrt(i_mpreal*i_mpreal) + cbrt((i_mpreal*i_mpreal*i_mpreal*i_mpreal)/(4*mju*mju)) );
			//sum += Factorial(i) / Factorial(k)  *  Factorial(i-2) / Factorial(i-k)  *  x.at(k) / Factorial(k-1)  *   x.at(i-k) / Factorial(i-k-1)  *  K;
		}
		
		mpreal xn = 0.5 * sum ;
		x.push_back(xn);
		
		cout << xn << endl << endl;
	}

	// Transforming sequence of {x_n} to {omega_n}
	for(int i=2; i<x.size(); i++)
		x.at(i) /= Factorial(i-1);

	cout << "{xn} done." << endl;

	
	// Now, I make sure that objects for calculating Bell polynomials are ready to use.

	for(int i=0; i<=M; ++i)
		for(int j=0; j<=M; ++j)
			done_Bells[i][j] = 0;

	for(int i=0; i<=M; ++i)
		for(int j=0; j<=M; ++j)
			calculated_Bells[i][j] = mpreal("0.0");

	for(int i=0; i<=M+1+1; ++i)
		for(int j=0; j<=M; ++j)
			done_coeffs[i][j] = 0;

	for(int i=0; i<=M+1+1; ++i)
		for(int j=0; j<=M; ++j)
			calculated_coeffs[i][j] = mpreal("0.0");


	 // Now, we start calculating the average number of clusters of given size, <n_g>.

	cout << "Calculating <n_g> starts." << endl;

	clock_t begin = clock();

	// I create vectors for <n_g> and for standard deviation.
	vector<mpreal> ng(max_g);
	vector<mpreal> sigma_g(max_g);

	// In this loop we take specific size g and calculate <n_g>.
	for(int g=1; g<=max_g; ++g)
	{
		cout << "g = " << g << endl;	

		// Expression for arbitrary kernel
		// In this expression, we used "x_g" but in fact it is "omega_g", as we have divided these values by (g-1)! above.
		ng.at(g-1) = BinomialCoefficient(M,g) * x.at(g) * Bell_g_k_specific_xn( M-g, k-1 ) / Bell_g_k_specific_xn( M, k );

		mpreal pomoc;

		// std dev
		if( 2*g > M ) pomoc = 0.;
		else
		{
			// <ns(ns-1)> for arbitrary kernel
			pomoc = BinomialCoefficient(M, g) * BinomialCoefficient(M-g, g) / Bell_g_k_specific_xn( M, k ) * x.at(g) * x.at(g) * Bell_g_k_specific_xn( M-2*g, k-2 ) ;
		}

		sigma_g.at(g-1) = mpfr::sqrt( pomoc + ng.at(g-1) - ng.at(g-1)*ng.at(g-1) );

	}

	cout << "Time elapsed: " << double(clock() - begin) / CLOCKS_PER_SEC << " s" << endl;

	// Saving results

	string filename_ng_str = "data/theory_const2_M=250_t=224_ng.txt";
	string filename_stddev_str = "data/theory_const2_M=250_t=224_stddev.txt";

	cout << "Saving <n_g> as " << filename_ng_str << endl;

	ofstream outputfile(filename_ng_str.c_str());
	for(auto n:ng)
		outputfile << n << endl;
	outputfile.close();

	cout << "Saving std. dev. as " << filename_stddev_str << endl;

	ofstream outputfile2(filename_stddev_str.c_str());
	for(auto s:sigma_g)
		outputfile2 << s << endl;
	outputfile2.close();


  return 0;
}

//********************** Factorial() ************************
// This function calculates factorial of a number.
mpreal Factorial(size_t n)
{
	mpreal factorial = 1;
	for (size_t i = 2; i <= n; i++)
	{
		factorial = factorial * i;
	}
	return factorial;
}

//******************** BinomialCoefficient() ****************
// This function calculates binomial coefficient.
mpreal BinomialCoefficient(int n, int k)
{
	return Factorial(n) / Factorial(k) / Factorial(n - k);
}

//******************* Bell_g_k_specific_xn() ****************
// This function calculates Bell polynomial with subscripts g and k.
// It uses {xn} sequence defined in the outer scope (lines 80-116).
mpreal Bell_g_k_specific_xn(int n, int k)
{
	int g = n - k + 1;

	//Computing the array of binominal coefficients

	vector< vector<mpreal> > binCoefArray;

	binCoefArray.resize(n + 1+1, vector<mpreal>(n + 1, 0));

	for (int i = 0; i <= n; i++)
	{
		binCoefArray[0][i] = 1;
	}
	for (int i = 0; i <= n; i++)
	{
		binCoefArray[1][i] = i;
	}

	int temp2;
	for (int i= 1; i <= n; i++)
	{
		if (g > i)
			temp2 = i;
		else if (i > g)
			temp2 = g - 1;
		else
			temp2 = i;

		for (int j = 1; j <= temp2 + 1; j++)
		{
			// check if coefficient for this Bell_g_k was calculated before
			if( done_coeffs[j][i] == 1 )
			{
				// if yes, then take it from memory and go on
				binCoefArray[j][i] = calculated_coeffs[j][i];
				continue;
			}

			// if it was not calculated before, calculate it and remember

			binCoefArray[j][i] = binCoefArray[j - 1][i] * (i - j + 1) / j;

			calculated_coeffs[j][i] = binCoefArray[j][i];
			done_coeffs[j][i] = 1;
		}
	}

	//Computing the Bell polynomial coefficients

	vector< vector<mpreal> > bell;
	bell.resize(k + 1, vector<mpreal>(n + 1, 0));
	bell[0][0] = 1;
	for (int j = 1; j <= k; j++)
	{
		for (int l = 1; l <= n; l++)
		{
			// check if this Bell_g_k was calculated before
			if( done_Bells[j][l] == 1 )
			{
				// if yes, then take it from memory and go on
				bell[j][l] = calculated_Bells[j][l];
				continue;
			}

			// if it was not calculated before, calculate it and remember

			for (int i = 1; i <= l - j + 1; i++)
			{
				bell[j][l] = bell[j][l] + binCoefArray[i - 1][l - 1] * bell[j - 1][l - i] * x[i];
			}

			calculated_Bells[j][l] = bell[j][l];
			done_Bells[j][l] = 1;
		}
	}

	//cout << endl << bell[k][n] << endl;

	return bell[k][n];
}
