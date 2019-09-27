#include <iostream>
#include <vector>
#include <string>
#include <fstream>
using namespace std;

#include "domain.hpp"
#include "datamesh.cpp"

TStepper::TStepper(DataMesh<double> &U, string method, string differencing, Patch &patch, int &time_steps, double &cf, double &cs, bool write_datafile)
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
      filename="data_files/ForwardEuler-"+differencing+"-SineWave-cf-"+to_string(cf)+".dat";
    }
    else if (method == "RungeKutta3")
    {
      filename="data_files/RungeKutta3-"+differencing+"-SineWave-cf-"+to_string(cf)+".dat";
    }
    datafile.open(filename);
    U.write(datafile);
  }

  DataMesh<double> dU = DataMesh<double>(dim,exts);
  for(int t=0;t<time_steps;t++)
  {
    ComputeRHS(U,dU,cs,dx,method,dt,differencing);
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
