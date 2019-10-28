#include <vector>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;



class RootFinder
{
  public:
    double Bisection(const double (*func)(const double &x), const double &a, 
                      const double &b, const double &tol, const int &maxIt);
    double Newton(const double (*func)(const double &x),
                  const double (*func_prime)(const double &x),
                  const double &x0, const double &tol, const int &maxIt);
    double Newton(const double (*func)(const double &x),
                  const double &epsilon, const double &x0, const double &tol,
                  const int &maxIt);
};
