#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <assert.h>

using namespace std;


vector<double> poly_approx_coeff (int n, double x)
{
	vector<double> poly_frac;
	for(int i=0;i<n;i++)
	{
		poly_frac.push_back(1.);
		for(int j=0;j<n;j++)
		{
			if (j!=i)
			{
				poly_frac[i] *= (x-j)/(i-j);
			}
		}
	}
	return poly_frac;
}

vector<double> poly_approx_coeff (int n, double x, double dx)
{
	vector<double> poly_frac;
	for(int i=0;i<n;i++)
	{
		poly_frac.push_back(1.);
		for(int j=0;j<n;j++)
		{
			if (j!=i)
			{
				poly_frac[i] *= (x-j*dx)/(i*dx-j*dx);
			}
		}
	}
	return poly_frac;
}

double poly_approx (int n, double x, vector<double> fi)
{
	double f_x = 0.;
	vector<double> poly_frac = poly_approx_coeff(n,x);
	for(int i=0;i<n;i++)
	{
		f_x += poly_frac[i]*fi[i];
	}
	return f_x;
}

double poly_approx (int n, double x, double dx, vector<double> fi)
{
	double f_x = 0.;
	vector<double> poly_frac = poly_approx_coeff(n,x,dx);
	for(int i=0;i<n;i++)
	{
		f_x += poly_frac[i]*fi[i];
	}
	return f_x;
}


double func_frac (double x)
{
	return (x*x-1.)/(x+2.);
}

double func_step(double x)
{
	return ( x >= 0.15 ? 1. : 0.);
}

int main()
{

	/* Part 2:
  	
	NOTE: Delta x (cell spacing)  := Dx
	      delta x (max(x,Dx))     := dx

	Functions defined above take into account two types of input,
	x = A*Dx (A=some constant), and x =/= A*Dx. If x is a scaled value of 
	dx, then Lagrange's formula simplifies. If not, then we must do the full calculations.
	
	To estimate f(2.5*Dx), assuming a given values of f_i, we run

	double f = poly_approx(int 5, double 2.5, vector<double> f_i);

	and to estimate f(5.5*Dx), we run

	double f = poly_approx(int 5, double 5.5, vector<double>,f_i);
	*/

	// ============================
	// Part 3.a:
	// ============================
	
{
	// Needed Constants
	double Dx = 0.1;
	int n = 5;
	int m = 200;

	// Values of known points and function
	vector<double> xi;
	vector<double> fi;
	xi.resize(n);
	fi.resize(n);
	for(int i=0;i<n;i++)
	{
		xi[i] = i*Dx;
		fi[i] = func_frac(i*Dx);
	}
	


	// Calculation of actual function 
	vector<double> x, f, f_lagr, dx, log_diff_f, log_dx;
	for(int i=0;i<m;i++)
	{
		x.push_back(i*n*Dx/m);
		f.push_back(func_frac(i*Dx));
		f_lagr.push_back(poly_approx(n,x[i],Dx,fi));
		dx.push_back(( x[i] >= Dx ? x[i] : Dx));
		log_diff_f.push_back(log(abs(f_lagr[i]-f[i])));
		log_dx.push_back(log(dx[i]));
	}	

	// Write Data to file
	ofstream diff_log_file;
	diff_log_file.open("3.a.gen.dat");
	for(int i=0;i<m;i++)
	{
		diff_log_file << log_dx[i] << "    " << log_diff_f[i] << "    " << f[i] << "    " << f_lagr[i] << endl;
	}
	diff_log_file.close();

	// ============================
	// Part 3.b:
	// ============================
	
	for (int i=0;i<n;i++)
	{
		fi[i] = func_step(i*Dx);
	}

	for(int i=0;i<m;i++)
	{
		f[i] = func_step(i*Dx);
		f_lagr[i] = poly_approx(n,x[i],Dx,fi);
		log_diff_f[i]=log(abs(f_lagr[i]-f[i]));
	}	
	
	ofstream diff_step_log;
	diff_step_log.open("3.b.gen.dat");
	for(int i=0;i<m;i++)
	{
		diff_step_log <<  log_dx[i] <<"    " << log_diff_f[i] << "    " << f[i] << "    " << f_lagr[i] << endl;
	}
	diff_step_log.close();

	system("gnuplot logplots.gp");

}
	
	// ============================
	// Part 4.a
	// ============================
{
	
	ofstream coeff;
	coeff.open("CoefficientCalc");
	for (int i=2;i<=10;i=i+2)
	{
		double x = ((i-1.)/2.);
		vector<double> ai = poly_approx_coeff(i,x);

		// Test case
		if (i==2)
		{
			assert (ai[0]==0.5 && ai[1]==0.5);
		}

		coeff << "For n=" << i << ", coefficients are: ";
		for (int j=0;j<i;j++)
		{
			coeff << ai[j] << ", ";
		}
		coeff << endl;
	}


	// ============================
	// Part 4.b
	// ============================
	
	for (int i=3;i<=11;i=i+2)
	{
		double x = i/2.;
		vector<double> ai = poly_approx_coeff(i,x);

		// Test Case
		if (i==3)
		{
			assert (ai[0]=-0.125 && ai[1]==0.75 && ai[2]==0.375); 
		}

		coeff << "For n=" << i << ", coefficients are: ";
		for (int j=0;j<i;j++)
		{
			coeff << ai[j] << ", ";
		}
		coeff << endl;
	}

	coeff.close();
	
}

}
