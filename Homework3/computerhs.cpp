#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "domain.hpp"
#include "datamesh.cpp"

double 
centered( double &xip,  double &xim)
{
  return (xip-xim)/2.;
}

double 
downstream( double &xip,  double &xi)
{
  return (xip-xi);
}

double 
upstream( double &xi,  double &xim)
{
  return (xi-xim);
}

double 
k1 ( double &xi,  double &xim, 
             double &cs,  double &dx,  double &dt)
{
  return xi-(dt/dx)*cs*downstream(xi,xim);
}

double 
k2 ( double &xi,  double &xim, 
             double &cs,  double &dx,  double &dt)
{
  return xi+k1(xi,xim,cs,dx,dt)+k1(xi,xim,cs,dx,dt);
}

ComputeRHS::ComputeRHS(DataMesh<double> &U, DataMesh<double> &dU,
                          double &cs, double &dx, string method, double &dt, string differencing)
{
  if(method=="ForwardEuler")
  {
    ForwardEuler(U,dU,cs,dx,dt,differencing);
  }
  else if (method=="RungeKutta3")
  {
    RungeKutta3(U,dU,cs,dx,dt);
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
          uim = U[i-1];
          val = U[i]-cs*(dt/dx)*centered(uip,uim);
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
          ui = U[i];
          val = U[i]-cs*(dt/dx)*downstream(uip,ui);
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
          ui = U[i];
          uim = U[i-1];
          val = U[i]-cs*(dt/dx)*upstream(ui,uim);
          dU.set_data_point(i,val);
        }
      }
    }
  }
}

void
ComputeRHS::RungeKutta3(const DataMesh<double> &U, DataMesh<double> &dU, const double &cs, const double &dx, const double &dt)
{
  k.resize(3);
  double ui,uim,val;

  for(int i=0;i<k.size();i++)
  {
    k[i].resize(U.get_total_points());

  }

  for(int i=0;i<k[0].size();i++)
  {

    if(i==0)
    {
      continue;
    }
    ui = U[i];
    uim = U[i-1];
    k[0][i]=U[i]+(dt/dx)*cs*downstream(ui,uim);

    if(i==1)
    {
      continue;
    }
    k[1][i] = U[i]+k[0][i]*(dt/2.)-(dt/dx)*cs
          *(downstream(ui,uim)+k[0][i]*(dt/2.)-k[0][i-1]*(dt/2.));
    if(i==2)
    {
      continue;
    }
    k[2][i] = U[i]+k[0][i]+2.*k[1][i]*dt-(dt/dx)*cs
          *(downstream(ui,uim)+2.*dt*(k[1][i]-k[1][i-1])-(k[0][i]-k[0][i-1]));

    if(! U.ghostzone(i))
    {
      val=U[i]+(1./6.)*dt*(k[0][i]+4.*k[1][i]+k[2][i]);
      dU.set_data_point(i,val);
    }
  }
}

