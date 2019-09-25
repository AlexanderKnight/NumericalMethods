#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "domain.hpp"
#include "datamesh.cpp"

ComputeRHS::ComputeRHS(Patch patch, DataMesh<double> &U, DataMesh<double> &dU,
                          double cs, string method, string differencing)
{
  if(method=="ForwardEuler")
  {
    ForwardEuler(patch,U,dU,cs,differencing);
  }
}

void
ComputeRHS::ForwardEuler(Patch patch, DataMesh<double> &U, DataMesh<double> &dU,
                          double cs, string differencing)
{
  double val;
  if(differencing=="centered")
  {
    for(int i=0;i<U.get_total_points();i++)
    {
      if (! U.ghostzone(i))
      {
        for(int dim=0;dim<U.get_dim();dim++)
        {
          val = -cs*(U.get_data_point(i+1)-U.get_data_point(i-1))/(2.*patch.dx(dim));
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
          val = -cs*(U.get_data_point(i+1)-U.get_data_point(i))/(patch.dx(dim));
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
          val = -cs*(U.get_data_point(i)-U.get_data_point(i-1))/(patch.dx(dim));
          dU.set_data_point(i,val);
        }
      }
    }
  }
}
