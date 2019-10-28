#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <cmath>
using namespace std;

#include "rf.hpp"


double
RootFinder::Bisection(const double (*func)(const double &x), const double &a,
                      const double &b, const double &tol, const int &maxIt)
{

  assert(a<b);
  assert(((*func)(a)*(*func)(b))<0.);

  double left = a;
  double right = b;

  double c;
  
  for (int i=0;i<maxIt;i++)
  {
    c = (left+right)/2.;
    if (abs((*func)(c))<tol||((right-left)/2.)<tol)
    {
      return c;
    }
    
    if((*func)(c)*(*func)(left)>0.)
    {
      left=c;
    }
    else
    {
      right=c;
    }
  }
  cout << "Failed to find root" << endl;
  return 0.;
}

double
RootFinder::Newton(const double (*func)(const double &x),
                    const double (*func_prime)(const double &x),
                    const double &x0, const double &tol, const int &maxIt)
{
  double xinit = x0;
  double x;

  for (int i=0;i<maxIt;i++)
  {
    x = xinit;
    xinit = x - ((*func)(x))/((*func_prime)(x));
    if((*func)(xinit)<tol)
    {
      return xinit;
    }
  }
  cout << "Failed to find root" << endl;
  return 0.;
}


double
RootFinder::Newton(const double (*func)(const double &x),
                   const double &epsilon, const double &x0, const double &tol,
                   const int &maxIt)
{
  double xinit = x0;
  double x;

  for(int i=0;i<maxIt;i++)
  {
    x=xinit;
    xinit = xinit-((*func)(x)*2.*epsilon)/((*func)(x+epsilon)-(*func)(x-epsilon));

    if((*func)(xinit)<tol)
    {
      return xinit;
    }
  }
  cout << "Failed to find root" << endl;
  return 0.;

}
