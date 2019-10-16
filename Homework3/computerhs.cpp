#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "domain.hpp"
#include "datamesh.cpp"


void
ComputeRHS::ForwardEuler(double (*fun)(const DataMesh<double> &, const int &, const double &),
                         const DataMesh<double> &U, DataMesh<double> &dU,
                         const double &cs, const double &dx, const double &dt)
  /* This currently only works for 1D problems. Soon to be reworked for 3D.
   */
{
  int Dim = U.get_dim();
  double val;
  for(int i=0;i<U.get_total_points();i++)
  {
    if (! U.ghostzone(i))
    {
      for(int dim=0;dim<Dim;dim++)
      {
        val = U[i]-cs*(dt/dx)*(*fun)(U,i,cs);
        dU.set_data_point(i,val);
      }
    }
  }
}



void
ComputeRHS::RungeKutta3(double (*fun)(const DataMesh<double> &, const int &, const double &cs), 
                        const DataMesh<double> &U, DataMesh<double> &dU, 
                        const double &cs, const double &dx, const double &dt)
{
  int dim = U.get_dim();
  vector<int> exts = U.get_exts();
  vector<int> gz_exts = U.get_gz_exts();
  bool periodic = true;
  double val;
  for(int i=0;i<exts.size();i++)
  {
    exts[i]-=2*gz_exts[i];
  }

  DataMesh<double> k1 = DataMesh<double>(dim,exts,gz_exts,periodic);
  DataMesh<double> k2 = DataMesh<double>(dim,exts,gz_exts,periodic);
  DataMesh<double> k3 = DataMesh<double>(dim,exts,gz_exts,periodic);
  DataMesh<double> k_temp = DataMesh<double>(dim,exts,gz_exts,periodic);
  for(int i=0;i<U.get_total_points();i++)
  {
    if(!k1.ghostzone(i))
    {
      val = (-1./dx)*((*fun)(U,i,cs));
      k1.set_data_point(i,val);
    }
  }
  k1.update_ghostzone();
 
  k_temp = k1;
  double dt_scale = dt/3.;
  k_temp *= dt_scale;
  k_temp += U;
  k_temp.update_ghostzone();
  for(int i=0;i<U.get_total_points();i++)
  {
    if(!k2.ghostzone(i))
    {
      val = (-1./dx)*((*fun)(k_temp,i,cs));
      k2.set_data_point(i,val);
    }
  }
  k2.update_ghostzone();

  // Evaluate k_temp = U + dt(-k1+2*k2)
  k_temp = k2;
  dt_scale = (2.*dt)/3.;
  k_temp *= dt_scale;
  k_temp += U;
  k_temp.update_ghostzone();
  for(int i=0;i<U.get_total_points();i++)
  {
    if(!k3.ghostzone(i))
    {
      val = (-1./dx)*((*fun)(k_temp,i,cs));
      k3.set_data_point(i,val);
    }
  }
  k3.update_ghostzone();
  //k3.print();
  for(int i=0;i<U.get_total_points();i++)
  {
    if(! U.ghostzone(i))
    {
      val=U[i]+(dt/4.)*(k1[i]+3.*k3[i]);
      dU.set_data_point(i,val);
    }
  }
}

