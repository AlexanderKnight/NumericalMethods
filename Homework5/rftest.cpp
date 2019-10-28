#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>
#include <cmath>
using namespace std;

#include "rf.hpp"

const double square(const double &x)
{
  return x*x-3.;
}

const double square_p(const double &x)
{
  return 2.*x;
}

const double cube(const double &x)
{
  return x*x*x-3.;
}

const double cube_p(const double &x)
{
  return 3.*x*x;
}

const double poly(const double &x)
{
  return 5.*pow(x,4)+12.*pow(x,3)-3.*pow(x,2)-6.*x+12;
}

const double poly_p(const double &x)
{
  return 20.*pow(x,3)+36.*pow(x,2)-6.*x-6.;
}

int main()
{
  double x0 = 4.;
  double a = 0.;
  double b = 3.;
  const double tol = 1.e-10;
  const int maxIt = 1000000;
  const int precis = 15;

  cout << endl <<"Parameters are:\nx0="<<x0<<"\na="<<a<<"\nb="<<b<<"\ntolerance="<<tol;
  cout << "\nMax Iterations="<<maxIt<<"\nPrecision="<<precis<<endl;
  cout << endl;

  cout << "Three test functions are:\nSquare: f(x)=-3+x^2\nCube: f(x)=-3+x^3";
  cout << "\nPolynomial: f(x)=5*x^4+12*x^3-3*x^2-6*x+12"<<endl;
  cout << endl;

  RootFinder rf;

  // Squares

  cout << "Square solutions: +/- sqrt(3) ~= +/- 1.7320508" << endl;
  double rfSqBi = rf.Bisection(square,a,b,tol,maxIt);
  assert(square(rfSqBi)<tol);
  cout << "Bisection root finder for square passes check"  << endl;
  cout << "Root=" << setprecision(precis) <<  rfSqBi << endl;

  double rfSqNe = rf.Newton(square,square_p,x0,tol,maxIt);
  assert(square(rfSqNe)<tol);
  cout << "Newton with known derivative root finder for square passes check"  << endl;
  cout << "Root=" <<  setprecision(precis) << rfSqNe << endl;

  double epsilon = 0.02;
  double rfSqNeDer = rf.Newton(square,epsilon,x0,tol,maxIt);
  assert(square(rfSqNeDer)<tol);
  cout << "Newton with unknown derivative root finder for square passes check" << endl;
  cout << "Root=" << setprecision(precis) <<  rfSqNeDer << endl << endl;

  // Cubes

  cout << "Cube solution: (3)^(1/3)~=1.44224957" << endl;
  double rfCuBi = rf.Bisection(cube,a,b,tol,maxIt);
  assert(cube(rfCuBi)<tol);
  cout << "Bisection root finder for cube passes check"  << endl;
  cout << "Root=" << setprecision(precis) <<  rfCuBi << endl;

  double rfCuNe = rf.Newton(cube,cube_p,x0,tol,maxIt);
  assert(cube(rfCuBi)<tol);
  cout << "Newton with known derivative root finder for cube passes check"  << endl;
  cout << "Root=" << setprecision(precis) <<  rfCuNe << endl;

  double rfCuNeDer = rf.Newton(cube,epsilon,x0,tol,maxIt);
  assert(cube(rfCuNeDer)<tol);
  cout << "Newton with unknown derivative root finder for cube passes check" << endl;
  cout << "Root=" << setprecision(precis) <<  rfCuNeDer << endl << endl;

  // Polynomials

  cout << "Polynomial solutions: -2.1990, -1.4449" << endl;
  a = -2.5;
  b = -1.5;
  epsilon=1.e-5;
  x0 = 0.;
  double rfPoBi = rf.Bisection(poly,a,b,tol,maxIt);
  assert(cube(rfPoBi)<tol);
  cout << "Bisection root finder for polynomial passes check"  << endl;
  cout << "Root=" << setprecision(precis) <<  rfPoBi << endl;

  double rfPoNe = rf.Newton(poly,poly_p,x0,tol,maxIt);
  assert(cube(rfPoBi)<tol);
  cout << "Newton with known derivative root finder for polynomial passes check"  << endl;
  cout << "Root=" << setprecision(precis) <<  rfPoNe << endl;

  double rfPoNeDer = rf.Newton(poly,epsilon,x0,tol,maxIt);
  assert(poly(rfPoNeDer)<tol);
  cout << "Newton with unknown derivative root finder for polynomial passes check" << endl;
  cout << "Root=" << setprecision(precis) <<  rfPoNeDer << endl << endl;
}


