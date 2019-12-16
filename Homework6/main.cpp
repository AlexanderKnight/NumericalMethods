#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <chrono>
#include <random>
using namespace std;

#include "domain.hpp"
#include "datamesh.cpp"
#include "rf.hpp"


int main()
{
  bool runPart1=1;
  bool runPart2=1;
  bool runPart3=1;
  bool runPart4=1;

  vector<vector<double>> test = {{1.,2.},{3.,4.}};
  vector<vector<double>> test_inv = inverse_matrix(test);
  cout << "Test Matrix: " << endl;
  for(int i=0;i<test.size();i++)
  {
    for(int j=0;j<test[i].size();j++)
    {
      cout << test[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
  cout << "Inverse Test Matrix: " << endl;
  for(int i=0;i<test_inv.size();i++)
  {
    for(int j=0;j<test_inv[i].size();j++)
    {
      cout << test_inv[i][j] << " ";
    }
    cout << endl;
  }

  // Part 1
  cout<<endl<<"----Part 1----"<<endl;

  cout.precision(13);
  double SqrtDetg = 1.;
  double rho0 = 1.e-4;
  double T_ref = 2.e-4;
  vector<double> u{0.7,0.0,0.};
  vector<vector<double>> g{{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
  double W_ref = 1.;
  for(int i=0;i<g.size();i++)
  {
    for(int j=0;j<g[0].size();j++)
    {
      W_ref += g[i][j]*u[i]*u[j];
    }
  }
  W_ref = sqrt(W_ref);
  
  double P_ref = 100.*rho0*rho0+rho0*T_ref;
  double rho_ref = SqrtDetg*W_ref*rho0;
  double h_ref = 1.+2.*P_ref/rho0;
  vector<double> Si (3,0.);
  for(int i=0;i<Si.size();i++)
  {
    Si[i] = rho_ref*h_ref*u[i];
  }
  double S_ref = 0.;
  for(int i=0;i<g.size();i++)
  {
    for(int j=0;j<g[0].size();j++)
    {
      S_ref += g[i][j]*Si[i]*Si[j];
    }
  }
  S_ref = sqrt(S_ref);
  double tau_ref = rho_ref*(h_ref*W_ref-1.)-SqrtDetg*P_ref;

  // Limits of W and T
  double WLimMin = 1.;
  double WLimMax = sqrt(1.+(S_ref/rho_ref)*(S_ref/rho_ref));
  double TLimMin = -0.1;
  double TLimMax = 1.;

  // Root-Finding variables
  double tol = 1.e-8;
  int maxIt = 1e4;

  //RootFinder rf=RootFinder(rho_ref, SqrtDetg,S_ref,tau_ref);
  RootFinder rf;
  if(runPart1)
  {
    auto start_t = chrono::high_resolution_clock::now();
    vector<double> guessT = rf.Find_Lorentz(WLimMin,WLimMax,TLimMin,TLimMax,
                                rho_ref,tau_ref,S_ref,SqrtDetg,tol,maxIt,'T');
    auto end_t = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsedT = end_t - start_t;
    cout << "Results from finding temp first are W=" << guessT[0] << ", T=" << guessT[1] << endl;
    cout << "Time for Temp run is "<< elapsedT.count() << " s" << endl;
  }

  // Part 2
  cout<<endl<<"----Part 2----"<<endl;

  if(runPart2)
  {
    auto start_p = chrono::high_resolution_clock::now();
    vector<double> guessP = rf.Find_Lorentz(WLimMin,WLimMax,TLimMin,TLimMax,
                                rho_ref,tau_ref,S_ref,SqrtDetg,tol,maxIt,'P');
    auto end_p = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsedP = end_p - start_p;

    cout << "Results from Pres are W=" << guessP[0] << ", T=" << guessP[1] << endl;
    cout << "Time for finding pressure first is "<< elapsedP.count() << " s" << endl;
  }
  
  // Part 3
  
  cout<<endl<<"----Part 3----" << endl;

  double T_guess = 1.9e-4; // correct value is 2.e-4
  double W_guess = 1.1;   // correct value is 1.0

  double epsilon = 1.e-12;

  if(runPart3)
  {
    auto start_N = chrono::high_resolution_clock::now();
    vector<double> guessTwoDNewton = rf.TwoDNewton(W_guess,T_guess,rho_ref,SqrtDetg,
                                              S_ref,tau_ref,epsilon,tol,maxIt);
    auto end_N = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsedN = end_N - start_N;
    cout << "Results from 2DNewton are W=" << guessTwoDNewton[0] << ", T=" << guessTwoDNewton[1] << endl;
    cout << "Time for 2DNewton run is "<< elapsedN.count() << " s" << endl;

    /* Note: While this is not currently faster than the previous methods,
     *        I think that this is due to the function to calculate the 
     *        determinant and inverse matrix of the Jacobian. I am going to
     *        leave it for the time being, since it is stable and calculating
     *        W and T correctly, but I will probably come back to fix this.
     */
  }

  // Part 4
  cout <<endl<< "----Part 4----" << endl;

  int testCount =100;
  //double testTol = 1.e-6;

  if(runPart4)
  {
    vector<bool> partsToRun = {1,1,1};
    vector<int> correct_counts = rf.RobustCheck(-12,-3,-6,-1,0,2,testCount,partsToRun,tol,maxIt,epsilon);
    cout << "The number of times that the root finding methods correctly determined the primative variables {Bisection-Temp, Bisection-Pressure, 2DNewton} : {";
    cout << correct_counts[0]<<", "<<correct_counts[1]<<", "<<correct_counts[2]<<"}"<<endl;
  }

}

