#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <chrono>
#include <random>
#include <assert.h>
using namespace std;

#include "domain.hpp"
#include "datamesh.cpp"
#include "rf.hpp"
#include "min.hpp"


// Functions needed for parts 1 & 2
const double cube (const double &x)
{
  return x*x*x;
} 
const double square(const double &x)
{
  return (x-2)*(x-2);
}
const double noise(const double &x)
{
  double alpha=-0.143*sin(1.75*(x+1.73));
  double beta=-0.18*sin(2.96*(x+4.98));
  double delta=0.012*sin(6.23*(x+3.17));
  double gamma=0.088*sin(8.07*(x+4.63));
  double val=alpha+beta+delta+gamma;
  return val;
}

// global variables for parts 3 & 4
// change value of p to change starting point of parts 3 & 4
// change value of xi to change ending point of part 3
extern vector<double> p = {-9.,-9.};
extern vector<double> xi = {9.,9.};
extern int n = p.size();
extern vector<double> xt(n);

// Functions needed for parts 3 & 4
const double well(const vector<double> &xii)
{
  return xii[0]*xii[0]+xii[1]*xii[1];
}

const double fpass(const double &x)
{
  for(int j=0;j<n;j++)
  {
    xt[j]=p[j]+x*(xi[j]-p[j]);
  }
  return well(xt);
}



int main()
{
  // Part 1 Test

  cout << setprecision(15);
  
  const double a = -100.;
  const double b = -99.;
  const double tol = 1.e-9;
  const int maxIt = 1e3;

  RootFinder rf = RootFinder();
  vector<double> bounds = rf.BracketMinimum(*square,a,b,tol,maxIt);
  cout<<"Bounds:"<<bounds[0]<<","<<bounds[1]<<","<<bounds[2]<<endl;
  assert(square(bounds[0])>square(bounds[1]) && square(bounds[1])<square(bounds[2]));

  // Part 2

  double min = rf.Brent(*square,bounds[0],bounds[1],bounds[2],tol,maxIt);
  cout << "Min of (x-2)^2 is " << min<<endl;

  const double noise_a = -1.;
  const double noise_b = -0.5;
  vector<double> noise_bounds = rf.BracketMinimum(*noise,noise_a,noise_b,tol,maxIt);
  cout<<"noise bounds: "<<noise_bounds[0]<<","<<noise_bounds[1]<<","<<noise_bounds[2]<<endl;
  double noise_min = rf.Brent(*noise,noise_bounds[0],noise_bounds[1],noise_bounds[2],tol,maxIt);
  cout<<"Min of noise is "<<noise_min<<endl;

  // Part 3

  double ax, xx, xmin;
  ax=0.0;
  xx=1.0;
  ofstream wellfile;
  wellfile.open("TwoDWellLineData.dat");
  double x_samp;
  for(int i=0;i<1000;i++)
  {
    x_samp = (double) i/(1000.);
    wellfile << x_samp << " " << fpass(x_samp) << "\n";
  }
  wellfile.close();

  vector<double> twoDBracket = rf.BracketMinimum(*fpass,ax,xx,tol,maxIt);
  cout<<"TwoDBracket is:"<<twoDBracket[0]<<", "<<twoDBracket[1]<<", "<<twoDBracket[2]<<endl;
  double twoDMin = rf.Brent(*fpass,twoDBracket[0],twoDBracket[1],twoDBracket[2],tol,maxIt);
  vector<double> xi_sol=xi;
  xi_sol[0]=xi[0]-p[0];
  xi_sol[1]=xi[1]-p[1];
  vector<double> p_sol=p;
  for(int j=0;j<n;j++)
  {
    xi_sol[j] *= twoDMin;
    p_sol[j] += xi_sol[j];
  }
  cout<<"Min of 2D plane line from ("<<p[0]<<","<<p[1]<<") to ("<<xi[0]<<","<<xi[1]<<") is ("<<p_sol[0]<<","<<p_sol[1]<<")"<<endl;

  // Part 4

  const double epsilon = 1e-12;
  const double tol2 = 1e-16;
  TwoDMinimum td = TwoDMinimum(*fpass,*well,*FRPR);
  vector<double> tdmin = td.FindMinimum(*fpass,*well,*FRPR,p,epsilon,tol2,maxIt);
  cout << "Min of 2D space is ("<<tdmin[0]<<","<<tdmin[1]<<")"<<endl;


}


