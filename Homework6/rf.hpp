#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
using namespace std;



class RootFinder
{
  public:
    double Bisection(const double (*func)(const double &x), const double &a, 
                      const double &b, const double &tol);
    double Newton(const double (*func)(const double &x),
                  const double (*func_prime)(const double &x),
                  const double &x0, const double &tol, const int &maxIt);
    double Newton(const double (*func)(const double &x),
                  const double &epsilon, const double &x0, const double &tol,
                  const int &maxIt);
    vector<double> Find_Lorentz(double &WLimMin, double &WLimMax, double &TLimMin, 
                  double &TLimMax, const double &rho, const double &tau,
                  const double &S, const double &SqrtDetg,
                  const double &tol, const int &maxIt, const char &solveVar);
    double Find_Pressure(const double &W);
    double Find_Temp(const double &W, double &mTLimMin, double &mTLimMax,
                      const double &tol, const int &maxIt);
                  
    double Lorentz(const double &W, const double &T);

    double Temp(const double &W, const double &T);
    double Solve_Temp(const double &P, const double &W);
    vector<double> TwoDNewton(const double &W0, const double &T0,
                              const double &mrho, const double &mSqrtDetg,
                              const double &mS, const double &mtau,
                              const double &epsilon, const double &tol,
                              const int &maxIt);
    vector<int> RobustCheck(const int rhoLimLow,const int rhoLimHigh,
                            const int TLimLow, const int TLimHigh,
                            const int WLimLow, const int WLimHigh,
                            const int &testCount,const vector<bool> &partsToRun,
                            const double &tol,const int &maxIt, 
                            const double &epsilon);

  private:
    double rho;
    double SqrtDetg;
    double S;
    double tau;
    double P(double mW, double mT){return 100.*rho0(mW,mT)*rho0(mW,mT)+rho0(mW,mT)*mT;}
    double rho0(double mW, double mT){return rho/(SqrtDetg*mW);}
    double h(double mW, double mT){return 1.+((2.*P(mW,mT))/rho0(mW,mT));}
};


inline double determinant(vector<vector<double>> A)
{
  assert(A.size()==A[0].size());
  if(A.size()==1)
  {
    return A[0][0];
  }
  else if(A.size()==2)
  {
    return A[0][0]*A[1][1]-A[1][0]*A[0][1];
  }
  else if(A.size()==3)
  {
    double a = A[0][0]*(A[0][1]*A[1][2]-A[1][2]*A[2][1]);
    double b = A[1][0]*(A[2][1]*A[0][2]-A[0][1]*A[2][2]);
    double c = A[2][0]*(A[0][1]*A[1][2]-A[1][1]*A[0][2]);
    return a+b+c;
  }
  return 0.;
}

inline vector<vector<double>> inverse_matrix(vector<vector<double>> A)
{

  assert(A.size()==A[0].size());
  double detA = determinant(A);
  /*if(A.size()==1)
  {
    vector<vector<double>> answer = {{1./A[0][0]}};
    return answer;
  }
  else if (A.size()==2)
  {*/
  vector<vector<double>> answer = {{A[1][1],-1.*A[0][1]},{-1.*A[1][0],A[0][0]}};
  
  for(int i=0;i<A.size();i++)
  {
    for(int j=0;j<A[i].size();j++)
    {
      answer[i][j] *= 1./detA;
    }
  }
  return answer;
  //}
}

