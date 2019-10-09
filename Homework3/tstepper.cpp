#include <iostream>
#include <vector>
#include <string>
#include <fstream>
using namespace std;

#include "domain.hpp"
#include "datamesh.cpp"

TStepper::TStepper(DataMesh<double> &U, string method, 
    double (*diff_func)(const DataMesh<double> &, const int &, const double &), 
    Patch &patch, int &time_steps, double &cf, double &cs, bool write_datafile)
{
  int dim = U.get_dim();
  vector<int> exts=U.get_exts();
  vector<int> gz_exts=U.get_gz_exts();
  bool is_periodic=true;

  double dx = patch.dx(0);
  double dt = cf*dx;
  U.update_ghostzone();


  ofstream datafile;
  if(write_datafile)
  {
    string filename;
    if(method == "ForwardEuler")
    {
      filename="data_files/ForwardEuler-SineWave-cf-"+to_string(cf)+"-"+to_string(exts[0]-2*gz_exts[0])+".dat";
    }
    else if (method == "RungeKutta3")
    {
      filename="data_files/RungeKutta3-SineWave-cf-"+to_string(cf)+"-"+to_string(exts[0]-2*gz_exts[0])+".dat";
    }
    datafile.open(filename);
    U.write(datafile);
  }

  ComputeRHS rhs;

  if(method=="ForwardEuler")
  {
    DataMesh<double> dU = DataMesh<double>(dim,exts);
    for(int t=0;t<time_steps;t++)
    {
      rhs.ForwardEuler(diff_func,U,dU,cs,dx,dt);
      U = dU;
      U.update_ghostzone();
      if(write_datafile)
      {
        U.write(datafile);
      }
    }
    if(write_datafile)
    {
      datafile.close();
    }
  }
  else if(method=="RungeKutta3")
  {
    DataMesh<double> dU = DataMesh<double>(dim,exts);
    for(int t=0;t<time_steps;t++)
    {
      rhs.RungeKutta3(diff_func,U,dU,cs,dx,dt);
      U = dU;
      U.update_ghostzone();
      if(write_datafile)
      {
        U.write(datafile);
      }
    }
    if(write_datafile)
    {
      datafile.close();
    }
  }
}
