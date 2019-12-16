#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <math.h>
using namespace std;

//#include "rf.hpp"

inline double
magnitude(const vector<double> A)
{
  double mag = 0.;
  for(int i=0;i<A.size();i++)
  {
    mag += A[i]*A[i];
  }
  mag = sqrt(mag);
  return mag;
}

/*
class F1dim
{
  private:
    const vector<double> &p;
    const vector<double> &xi;
    int n;
    const double (*func)(const vector<double> &x);
    vector<double> xt;
  public:
    F1dim(vector<double> &pp, vector<double> &xii, 
          const double (*funcc)(const vector<double> &x));
    const static double f(const double &x);
};
*/

class TwoDMinimum
{
  private:
    const double (*func)(const double &x);
    const double (*TwoDFunc)(const vector<double> &x);
    const double (*beta)(const vector<double> &xni, const vector<double> &xnj);
    void steepest_descent(const double (*funcc)(const vector<double> &x),
                                    const vector<double> &x0, vector<double> &delx,
                                    const double epsilon);
  public:
    TwoDMinimum(const double (*func)(const double &x),
                const double (*TwoDFunc)(const vector<double> &x),
                const double (*beta)(const vector<double> &xni,
                                const vector<double> &xnj));
  
  const vector<double> FindMinimum(const double (*func)(const double &x),
                              const double (*TwoDFunc)(const vector<double> &x),
                              const double (*beta)(const vector<double> &xni,
                                                    const vector<double> &xnj),
                              const vector<double> &pp, const double &epsilon,
                              const double &tol, const int &maxIt);
};



inline const double
FRPR (const vector<double> &xni, const vector<double> &xnj)
{
  assert(xni.size()==xnj.size());
  double num=0., denom=0.;
  for(int i=0;i<xni.size();i++)
  {
    num += xnj[i]*(xnj[i]-xni[i]);
    denom += xni[i]*xni[i];
  }
  return max(0.,num/denom);
}
  
