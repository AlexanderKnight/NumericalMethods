#include <iostream>
#include <vector>
#include <string>
#include <fstream>
using namespace std;

#include "domain.hpp"
#include "datamesh.cpp"

TStepper::TStepper(DataMesh<double> &U, string method, 
    const double (*spacial_deriv)(const DataMesh<double> &, const int &, const double &, const double &), 
    const double (*wave_eq)(const double &, const double &, const double &),
    double &dx, int &time_steps, double &cf, double &cs, bool write_datafile)
{
  int dim = U.get_dim();
  vector<int> exts=U.get_exts();
  vector<int> gz_exts=U.get_gz_exts();
  bool is_periodic=true;

  double dt = cf*dx;
  U.update_ghostzone();


  // Sets up for writing data to file
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
    else if (method == "LLF")
    {
      filename="data_files/LLF-SineWave-cf-"+to_string(cf)+"-"+to_string(exts[0]-2*gz_exts[0])+".dat";
    }
    datafile.open(filename);
    U.write(datafile);
  }

  ComputeRHS rhs;
  DataMesh<double> dU = DataMesh<double>(dim,exts);

  if(method=="ForwardEuler")
  {
    for(int t=0;t<time_steps;t++)
    {
      rhs.ForwardEuler(spacial_deriv,wave_eq,U,dU,cs,dx,dt);
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
    for(int t=0;t<time_steps;t++)
    {
      U.update_ghostzone();
      rhs.RungeKutta3(spacial_deriv,wave_eq,U,dU,cs,dx,dt);
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
  else if(method=="LLF")
  {
    for(int t=0;t<time_steps;t++)
    {
      U.update_ghostzone();
      rhs.LLF(wave_eq,U,dU,cs,dx,dt);
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
