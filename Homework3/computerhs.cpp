#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "domain.hpp"
#include "datamesh.cpp"

double 
centered(double &xip,  double &xi, double &xim)
{
  return (xip-xim)/2.;
}

double 
downstream(double &xip,  double &xi, double &xim)
{
  return (xip-xi);
}

double 
upstream(double &xip, double &xi,  double &xim)
{
  return (xi-xim);
}

/*
double 
k1 (double &xip, double &xi,  double &xim, 
             double &cs,  double &dx,  double &dt)
{
  return xi-(dt/dx)*cs*downstream(uip,xi,xim);
}

double 
k2 (double &xip, double &xi,  double &xim, 
             double &cs,  double &dx,  double &dt)
{
  return xi+k1(xip,xi,xim,cs,dx,dt)+k1(xip,xi,xim,cs,dx,dt);
}
*/

ComputeRHS::ComputeRHS(DataMesh<double> &U, DataMesh<double> &dU,
                          double &cs, double &dx, string method, double &dt, string differencing)
{
  if(method=="ForwardEuler")
  {
    ForwardEuler(U,dU,cs,dx,dt,differencing);
  }
  else if (method=="RungeKutta3")
  {
    RungeKutta3(U,dU,differencing,cs,dx,dt);
  }
}

void
ComputeRHS::ForwardEuler(DataMesh<double> &U, DataMesh<double> &dU,
                          double &cs, double &dx, double &dt,string differencing)
  /* This currently only works for 1D problems. Soon to be reworked for 3D.
   */
{
  double val,uip,ui,uim;
  if(differencing=="centered")
  {
    for(int i=0;i<U.get_total_points();i++)
    {
      if (! U.ghostzone(i))
      {
        for(int dim=0;dim<U.get_dim();dim++)
        {
          uip = U[i+1];
          ui  = U[i];
          uim = U[i-1];
          val = U[i]-cs*(dt/dx)*centered(uip,ui,uim);
          dU.set_data_point(i,val);
        }
      }
    }
  }
  if(differencing=="downstream")
  {
    for(int i=0;i<U.get_total_points();i++)
    {
      if (! U.ghostzone(i))
      {
        for(int dim=0;dim<U.get_dim();dim++)
        {
          uip = U[i+1];
          ui  = U[i];
          uim = U[i-1];
          val = U[i]-cs*(dt/dx)*downstream(uip,ui,uim);
          dU.set_data_point(i,val);
        }
      }
    }
  }
  if(differencing=="upstream")
  {
    for(int i=0;i<U.get_total_points();i++)
    {
      if (! U.ghostzone(i))
      {
        for(int dim=0;dim<U.get_dim();dim++)
        {
          uip = U[i+1];
          ui  = U[i];
          uim = U[i-1];
          val = U[i]-cs*(dt/dx)*upstream(uip,ui,uim);
          dU.set_data_point(i,val);
        }
      }
    }
  }
}

void
ComputeRHS::RungeKutta3(const DataMesh<double> &U, DataMesh<double> &dU, string diff,
                        const double &cs, const double &dx, const double &dt)
{
  k.resize(3);
  double uip,ui,uim,val;

  for(int i=0;i<k.size();i++)
  {
    k[i].resize(U.get_total_points());

  }

  for(int i=0;i<k[0].size()-1;i++)
  {

    if(i==0)
    {
      continue;
    }
    uip = U[i+1];
    ui = U[i];
    uim = U[i-1];
    if(diff=="downstream")
    {
      k[0][i]=(dt/dx)*cs*downstream(uip,ui,uim);
    }
    else if (diff=="upstream")
    {
      k[0][i]=(dt/dx)*cs*upstream(uip,ui,uim);
    }
    else if (diff=="centered")
    {
      k[0][i]=(dt/dx)*cs*centered(uip,ui,uim);
    }


    if(i==1)
    {
      continue;
    }
    if(diff=="downstream")
    {
    k[1][i] = k[0][i]*(dt/2.)-(dt/dx)*cs
          *(downstream(uip,ui,uim)+k[0][i+1]*(dt/2.)-k[0][i]*(dt/2.));
    }
    else if (diff=="upstream")
    {
    k[1][i] = k[0][i]*(dt/2.)-(dt/dx)*cs
          *(upstream(uip,ui,uim)+k[0][i]*(dt/2.)-k[0][i-1]*(dt/2.));
    }
    else if (diff=="centered")
    {
    k[1][i] = k[0][i]*(dt/2.)-(dt/dx)*cs
          *(centered(uip,ui,uim)+k[0][i+1]*(dt/2.)-k[0][i-1]*(dt/2.));
    }


    if(i==2)
    {
      continue;
    }
    if(diff=="downstream")
    {
    k[2][i] = k[0][i]+2.*k[1][i]*dt-(dt/dx)*cs
          *(downstream(uip,ui,uim)+2.*dt*(k[1][i+1]-k[1][i])-(k[0][i+1]-k[0][i]));
    }
    else if (diff=="upstream")
    {
    k[2][i] = k[0][i]+2.*k[1][i]*dt-(dt/dx)*cs
          *(upstream(uip,ui,uim)+2.*dt*(k[1][i]-k[1][i-1])-(k[0][i]-k[0][i-1]));
    }
    else if (diff=="centered")
    {
    k[2][i] = k[0][i]+2.*k[1][i]*dt-(dt/dx)*cs
          *(centered(uip,ui,uim)+2.*dt*(k[1][i+1]-k[1][i-1])-(k[0][i+1]-k[0][i-1]));
    }

    if(! U.ghostzone(i))
    {
      val=U[i]+(1./6.)*dt*(k[0][i]+4.*k[1][i]+k[2][i]);
      dU.set_data_point(i,val);
    }
  }
}

