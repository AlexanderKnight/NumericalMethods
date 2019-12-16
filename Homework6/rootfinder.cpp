#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <cmath>
#include <random>
using namespace std;

#include "rf.hpp"

double get_rand_double(int lowPow, int highPow)
{
  random_device rd;
  mt19937 rng(rd());
  uniform_real_distribution<double> uniDoub(0.00,9.99);
  uniform_int_distribution<int> uniInt(lowPow,highPow);
  auto rand_double = uniDoub(rng);
  auto rand_power_int = uniInt(rng);
  double rand_pow = pow(10.,rand_power_int);
  double rand = rand_double*rand_pow;
  return rand;
}



double
RootFinder::Temp(const double &mW, const double &mT)
{
  //cout<<"rf.Temp[W,T]: "<<W<<", "<<T<<endl;
  double answer = h(mW,mT)*mW-1.-(SqrtDetg*P(mW,mT)/rho)-(tau/rho);
  return answer;
}

// Might need work, double check
double
RootFinder::Solve_Temp(const double &P_temp, const double &mW)
{
  mW;
  double answer = ((P_temp-100.*(rho*rho)/(SqrtDetg*SqrtDetg*mW*mW))/(rho/(SqrtDetg*mW)));
  return answer;
}

double
RootFinder::Lorentz(const double &mW, const double &mT)
{
  //cout<<"rf.Lorentz[W,T]: "<<W<<", "<<T<<endl;
  double answer = mW*mW-1.-((S*S)/(rho*rho*h(mW,mT)*h(mW,mT)));
  return answer;
}

double
RootFinder::Find_Temp(const double &W,
                      double &mTLimMin, double &mTLimMax,
                      const double &tol, const int &maxIt)
{
  assert(mTLimMin<mTLimMax);
  //assert(TempMin*TempMax < 0.);

  double TLimMin = mTLimMin;
  double TLimMax = mTLimMax;

  double TLimMid, TempMin, TempMid, TempMax;
  for(int i=0;i<maxIt;i++)
  {
    TLimMid = (TLimMin+TLimMax)/2.;
    TempMid = Temp(W,TLimMid);
    TempMin = Temp(W,TLimMin);
    TempMax = Temp(W,TLimMax);
    if (abs((TLimMax-TLimMin)/2.)<tol)
    {
      return TLimMid;
    }

    if(TempMin*TempMax>0.) 
    {
      return -1;
    }
    else if(TempMin*TempMid > 0.)
    {
      TLimMin = TLimMid;
    }
    else
    {
      TLimMax = TLimMid;
    }
  }
  return -1.;
}

double
RootFinder::Find_Pressure(const double &W)
/*
 * Gives pressure for a given W, solves
 * h*W-1-(sqrt(gamma)*P/rho)-(tau/rho)=0
 * for P, with h=1+2P/rho0
 */
{
  double answer = ((rho-W*rho+tau)/(2.*W*W-SqrtDetg));
  return answer;
}

vector<double>
RootFinder::Find_Lorentz(double &mWLimMin, double &mWLimMax, double &mTLimMin, 
                          double &mTLimMax, const double &mrho, const double &mtau,
                          const double &mS, const double &mSqrtDetg,
                          const double &tol, const int &maxIt, const char &solveVar)
{
  rho=mrho;
  SqrtDetg=mSqrtDetg;
  S=mS;
  tau=mtau;
  assert(mWLimMin<mWLimMax);
  double WLimMin, WLimMax, WLimMid, WMin, WMax, WMid;
  //double TLimMin, TLimMax, WLimMin, WLimMax;
  WLimMin = mWLimMin;
  WLimMax = mWLimMax;

  
  //double WLimMid, WMin, WMax, WMid;
  //double WLimMid, TMin, TMax, TMid, WMin, WMax, WMid;
  vector<double> error = {-1.,-1.};



  if(solveVar=='T')
  {
    double TLimMin, TLimMax, TMin, TMax, TMid;
    for(int i=0;i<maxIt;i++)
    {
      WLimMid = (WLimMin+WLimMax)/2.; 
      TMin = Find_Temp(WLimMin,mTLimMin,mTLimMax,tol,maxIt);
      if(TMin==-1.) return error;
      TMax = Find_Temp(WLimMax,mTLimMin,mTLimMax,tol,maxIt);
      if(TMax==-1.) return error;
      TMid = Find_Temp(WLimMid,mTLimMin,mTLimMax,tol,maxIt);
      if(TMid==-1.) return error;

      WMin = Lorentz(WLimMin,TMin);
      if(WMin==-1.) return error;
      WMax = Lorentz(WLimMax,TMax);
      if(WMax==-1.) return error;
      WMid = Lorentz(WLimMid,TMid);
      if(WMid==-1.) return error;

      
      if((abs((WLimMax-WLimMin)/2.)<tol) && ((abs(TMax-TMin)/2.)<tol))
      {
        vector<double> result = {WLimMid, TMid};
        //result.push_back(WLimMid);
        //result.push_back(TMid);
        return result;
      }

      if(WMin*WMax>0.)
      {
        //vector<double> error{-1.,-1.};
        return error;
      }
      else if(WMin*WMid>0.)
      {
        WLimMin = WLimMid;
      }
      else
      {
        WLimMax = WLimMid;
      }
    }
  }
  else if(solveVar=='P')
  {
    double PMin, PMax, PMid, TMin, TMax, TMid;
    for(int i=0;i<maxIt;i++)
    {
      WLimMid = (WLimMin+WLimMax)/2.; 
      PMin = Find_Pressure(WLimMin);
      PMax = Find_Pressure(WLimMax);
      PMid = Find_Pressure(WLimMid);

      TMin = Solve_Temp(PMin,WLimMin);
      if(TMin==-1.) return error;
      TMid = Solve_Temp(PMid,WLimMid);
      if(TMid==-1.) return error;
      TMax = Solve_Temp(PMax,WLimMax);
      if(TMax==-1.) return error;


      WMin = Lorentz(WLimMin,TMin);
      if(WMin==-1.) return error;
      WMax = Lorentz(WLimMax,TMax);
      if(WMax==-1.) return error;
      WMid = Lorentz(WLimMid,TMid);
      if(WMid==-1.) return error;

      
      
      if(abs((WLimMax-WLimMin)/2.)<tol)
      {
        vector<double> result;
        result.push_back(WLimMid);
        result.push_back(TMid);
        return result;
      }

      if(WMin*WMid>0.&&WMid*WMax>0.)
      {
        vector<double> error{-1.,-1.};
        return error;
      }

      else if(WMin*WMid>0.)
      {
        WLimMin = WLimMid;
      }
      else
      {
        WLimMax = WLimMid;
      }
    }
  }

  // Returns error if nothing found
  //vector<double> error{-1.,-1.};
  return error;
}

double
RootFinder::Bisection(const double (*func)(const double &x), const double &a,
                      const double &b, const double &tol)
{

  assert(a<b);
  assert(((*func)(a)*(*func)(b))<0.);

  double left = a;
  double right = b;

  double c;
  
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

vector<double>
RootFinder::TwoDNewton(const double &W0, const double &T0,
                       const double &mrho, const double &mSqrtDetg,
                       const double &mS, const double &mtau,
                       const double &epsilon, const double &tol,
                       const int &maxIt)
{
  rho=mrho;
  SqrtDetg=mSqrtDetg;
  S=mS;
  tau=mtau;
  //double W=W0;
  //double T=T0;
  
  double W_sol=W0;
  double T_sol=T0;
  double W_epp, W_epm, T_epp, T_epm;
  vector<double> F(2);
  vector<vector<double>> J(2,vector<double>(2));
  vector<vector<double>> InvJ(2,vector<double>(2));
  for(int i=0;i<maxIt;i++)
  {
    //W=W_sol;
    //T=T_sol;
    F[0] = Lorentz(W_sol,T_sol);
    F[1] = Temp(W_sol,T_sol);


    if(abs(F[0])<tol and abs(F[1])<tol)
    {
      vector<double> answer = {W_sol,T_sol};
      return answer;
    }
      

    // Taking the Jacobian to be
    //
    // J =| d_w(Lorentz(w,t))   d_t(Lorentz(w,t))| 
    //    | d_w(Temp(w,t))      d_t(Temp(w,t))   |
    //

    W_epp = W_sol+epsilon;
    W_epm = W_sol-epsilon;
    T_epp = T_sol+epsilon;
    T_epm = T_sol-epsilon;
    J[0][0] = (Lorentz(W_epp,T_sol)
                    -Lorentz(W_epm,T_sol))/(2.*epsilon);
    J[1][0] = (Temp(W_epp,T_sol)
                    -Temp(W_epm,T_sol))/(2.*epsilon);
    J[0][1] = (Lorentz(W_sol,T_epp)
                    -Lorentz(W_sol,T_epm))/(2.*epsilon);
    J[1][1] = (Temp(W_sol,T_epp)
                    -Temp(W_sol,T_epm))/(2.*epsilon);

    InvJ = inverse_matrix(J);

    W_sol -= InvJ[0][0]*F[0]+InvJ[0][1]*F[1];
    T_sol -= InvJ[1][0]*F[0]+InvJ[1][1]*F[1];
  }
  vector<double> answer = {-1.,-1.};
  return answer;
}

vector<int>
RootFinder::RobustCheck(const int rhoLimLow, const int rhoLimHigh,
                        const int TLimLow, const int TLimHigh,
                        const int WLimLow, const int WLimHigh,
                        const int &testCount,const vector<bool> &partsToRun,
                        const double &tol,const int &maxIt, const double &epsilon)
{
  cout.precision(13);
  //double SqrtDetg = 1.;
  // Recalculate this for the new W, see paper, talked to Sasha
  //vector<double> u{0.,0.,0.};
  vector<double> u{0.,0.,0.};
  vector<vector<double>> g{{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};

  vector<int> correct_counts = {0,0,0};

  // Limits of W and T
  double WLimMin_rand = 1.;
  //double WLimMax_rand = sqrt(1.+(S_rand_ref/rho_rand_ref)*(S_rand_ref/rho_rand_ref));
  double WLimMax_rand = 1000.;
  double TLimMin_rand = -1.e5;
  double TLimMax_rand = 1.e5;

  for(int count=0;count<testCount;count++)
  {
    // Get Random values for rho, T, W
    double rho0_rand = get_rand_double(rhoLimLow,rhoLimHigh);
    double T_rand = get_rand_double(TLimLow,TLimHigh);
    double W_rand = get_rand_double(WLimLow,WLimHigh)+1.;
    double rho_rand = rho0_rand*SqrtDetg*W_rand;


    u[0]=sqrt(W_rand*W_rand-1);

    // Recalculate the paramters needed for the algorithms
    double P_rand = 100.*rho0_rand*rho0_rand+rho0_rand*T_rand;
    //double rho_rand_ref = SqrtDetg*W_rand*rho_rand;
    double h_rand = 1.+2.*P_rand/rho0_rand;
    vector<double> Si_rand (3,rho_rand*h_rand);
    for(int i=0;i<Si_rand.size();i++)
    {
      Si_rand[i] *= u[i];
    }
    double S_rand_ref = 0.;
    for(int i=0;i<g.size();i++)
    {
      for(int j=0;j<g[0].size();j++)
      {
        S_rand_ref += g[i][j]*Si_rand[i]*Si_rand[j];
      }
    }
    S_rand_ref = sqrt(S_rand_ref);
    double tau_rand_ref = rho_rand*(h_rand*W_rand-1.)-SqrtDetg*P_rand;


    rho=rho_rand;
    SqrtDetg=1.;
    S=S_rand_ref;
    tau=tau_rand_ref;
    


    if(partsToRun[0])
    {
      vector<double> guessT_rand = Find_Lorentz(WLimMin_rand,WLimMax_rand,TLimMin_rand,TLimMax_rand,
                                  rho_rand,tau_rand_ref,S_rand_ref,SqrtDetg,tol,maxIt,'T');
      //cout <<"LorentzTDiff[W,T]: ["<<abs(guessT_rand[0]-W_rand)<<", "<<abs(guessT_rand[1]-T_rand)<<"]"<<endl;
     /* if(guessT_rand[0]==-1.||guessT_rand[1]==-1.)
      {
        continue;
      }*/

      if((abs(guessT_rand[0]-W_rand)<(1.e4*tol)) and (abs(guessT_rand[1]-T_rand)<1.e4*tol))
      {
        correct_counts[0]++;
      }
    }
    if(partsToRun[1])
    {
      vector<double> guessP_rand = Find_Lorentz(WLimMin_rand,WLimMax_rand,TLimMin_rand,TLimMax_rand,
                                  rho_rand,tau_rand_ref,S_rand_ref,SqrtDetg,tol,maxIt,'P');
      //cout <<"LorentzPDiff[W,T]: ["<<abs(guessP_rand[0]-W_rand)<<", "<<abs(guessP_rand[1]-T_rand)<<"]"<<endl;
      //if(guessP_rand[0]==-1.||guessP_rand[1]==-1.)
      //{
       // continue;
      //}
      if((abs(guessP_rand[0]-W_rand)<1.e4*tol) and (abs(guessP_rand[1]-T_rand)<1.e4*tol))
      {
        correct_counts[1]++;
      }
    }

    

    double T_rand_guess = 0.9*T_rand; 
    double W_rand_guess = 1.1*W_rand;   

    if(partsToRun[2])
    {
      vector<double> guessTwoDNewton_rand = TwoDNewton(W_rand_guess,T_rand_guess,rho_rand,SqrtDetg,
                                                S_rand_ref,tau_rand_ref,epsilon,tol,maxIt);
      //cout <<"TwoDNewtonDiff[W,T]: ["<<abs(guessTwoDNewton_rand[0]-W_rand)<<", "<<abs(guessTwoDNewton_rand[1]-T_rand)<<"]"<<endl;
      if((abs(guessTwoDNewton_rand[0]-W_rand)<1.e4*tol) and (abs(guessTwoDNewton_rand[1]-T_rand)<1.e4*tol))
      {
        correct_counts[2]++;
      }
    }
  }
  return correct_counts;
}


