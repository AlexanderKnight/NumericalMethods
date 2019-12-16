#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <cmath>
#include <random>
using namespace std;

#include "min.hpp"
#include "rf.hpp"

/*
F1dim::F1dim(vector<double> &pp, vector<double> &xii,
             const double (*funcc)(const vector<double> &x))
{
  p=pp;
  xi=xii;
  func=funcc;
  n=p.size();
  xt.resize(n);
}
*/

TwoDMinimum::TwoDMinimum(const double (*funcc)(const double &x),
                         const double (*TwoDFuncc)(const vector<double> &x),
                          const double (*betaa)(const vector<double> &xni,
                                                const vector<double> &xnj))
{
  func=funcc;
  TwoDFunc=TwoDFuncc;
  beta=betaa;
}

const vector<double>
TwoDMinimum::FindMinimum(const double (*func)(const double &xx),
                         const double (*TwoDFunc)(const vector<double> &x),
                          const double (*beta)(const vector<double> &xni,
                                          const vector<double> &xnj),
                          const vector<double> &pp, const double &epsilon,
                          const double &tol, const int &maxIt)
{
  vector<double> delx0 = pp;
  vector<double> delx1 = pp;
  const int n = pp.size();
  vector<double> xhold(n);
  vector<double> xi(n);


  steepest_descent(*TwoDFunc,pp,xi,epsilon);
  RootFinder rf = RootFinder();
  double ax=0.,xx=1.;
  vector<double> bracket = rf.BracketMinimum(*func,ax,xx,tol,maxIt);
  double min = rf.Brent(*func,bracket[0],bracket[1],bracket[2],tol,maxIt);

  vector<double> xi_sol=xi;
  for(int j=0;j<n;j++)
  {
    xi_sol[j] *= min;
    delx1[j] += xi_sol[j];
  }
  vector<double> s0 = xi;

  double bet;
  for(int iter=0;iter<maxIt;iter++)
  {
    steepest_descent(*TwoDFunc,delx1,xi,epsilon);
    bet=beta(delx0,delx1);
    s0[0] = xi[0]+bet*s0[0];
    s0[1] = xi[1]+bet*s0[1];
    xi = s0;
    bracket=rf.BracketMinimum(*func,ax,xx,tol,maxIt);
    min=rf.Brent(*func,bracket[0],bracket[1],bracket[2],tol,maxIt);

    delx0 = delx1;
    xi_sol=xi;
    xi_sol[0]-=delx1[0];
    xi_sol[1]-=delx1[1];
    for(int j=0;j<n;j++)
    {
      xi_sol[j] *= min;
      delx1[j] += xi_sol[j];
    }
    if(fabs(magnitude(delx0)-magnitude(delx1))<tol)
    {
      return delx1;
    }
  }
  cout << "No minimum found" << endl;
}


void
TwoDMinimum::steepest_descent(const double (*funcc)(const vector<double> &x),
                              const vector<double> &x0, vector<double> &delx,
                              const double epsilon)
{
  const int dimm = x0.size();
  vector<double> x0_temp(dimm);
  double ftempm;
  double ftempp;

  for(int dim=0;dim<dimm;dim++)
  {
    x0_temp = x0;
    x0_temp[dim] -= epsilon;
    ftempm=funcc(x0_temp);
    x0_temp = x0;
    x0_temp[dim] += epsilon;
    ftempp=funcc(x0_temp);
    delx[dim] = -(ftempp-ftempm)/(2.*epsilon);
  }
}

