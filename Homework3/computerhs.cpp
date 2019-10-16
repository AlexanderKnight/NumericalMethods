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
  /* This currently only works for 1D problems. Future work will change this
   * to 3D compatible.
   *
   * Note: (*fun): spatial derivative function, defined in domain.hpp
   */
{
  // Gets data to make k1, k2, k3,and k_temp for the calculations
  int dim = U.get_dim();
  vector<int> exts = U.get_exts();
  vector<int> gz_exts = U.get_gz_exts();
  bool periodic = true;
  double val;

  // By default the extents of a DataMesh include the ghostzones, but
  // the original extents are needed, so this removes the ghostzones
  // from the count.
  for(int i=0;i<exts.size();i++)
  {
    exts[i]-=2*gz_exts[i];
  }

  DataMesh<double> k1 = DataMesh<double>(dim,exts,gz_exts,periodic);
  DataMesh<double> k2 = DataMesh<double>(dim,exts,gz_exts,periodic);
  DataMesh<double> k3 = DataMesh<double>(dim,exts,gz_exts,periodic);
  DataMesh<double> k_temp = DataMesh<double>(dim,exts,gz_exts,periodic);

  /* Current method utilizes Kutta's Third order method with
   *  0  |  0   0   0
   * 1/2 | 1/2  0   0
   *  1  | -1   2   0
   * -----------------
   *     | 1/6 2/3 1/6
   *
   * Or that k1 = (d/dx)(U_x), k2 = (d/dx)(U_x+(dt/2)*k1), 
   * and k3=(d/dx)(U_x+2*k2-k1)
   */

  // Assign k1 t0 (d/dx)(f(x))
  for(int i=0;i<U.get_total_points();i++)
  {
    if(!k1.ghostzone(i))
    {
      val = (-1./dx)*((*fun)(U,i,cs));
      k1.set_data_point(i,val);
    }
  }
  k1.update_ghostzone();
 
  // Evaluate k_temp = U + (dt/2)*k1
  k_temp = k1;
  double dt_scale = dt/2.;
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
  k_temp = k1;
  k_temp *= -1.;
  k_temp += k2;
  k_temp += k2;
  k_temp *= dt;
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


  for(int i=0;i<U.get_total_points();i++)
  {
    if(! U.ghostzone(i))
    {
      val=U[i]+(dt/6.)*(k1[i]+4.*k2[i]+k3[i]);
      dU.set_data_point(i,val);
    }
  }
}

