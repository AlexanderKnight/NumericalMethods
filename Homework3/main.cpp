#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;

#include "domain.hpp"
#include "datamesh.cpp"

int main()
{
  int dimensions = 1;
  vector<int> extents = {1000};
  vector<int> gz_extents = {2};
  vector<vector<double> > limits = {{0.,2.*M_PI}};
  bool is_periodic=true;
  double cs=1.;
  int time_steps =10000;
  bool write_datafile = true;
  
  vector<double> sinx;
  sinx.resize(extents[0]);
  for(int i=0;i<sinx.size();i++)
  {
    sinx[i]=sin(i*(limits[0][1]-limits[0][0])/(static_cast<double>(extents[0])));
  }
  DataMesh<double> U = DataMesh<double>(dimensions, extents,gz_extents,is_periodic);
  Patch P = Patch(dimensions, extents, limits);
  double dx = P.dx(0);
  //vector<string> Method = {"ForwardEuler", "RungeKutta3"};
  //vector<string> differencing = {"centered", "downstream", "upstream"};
  //vector<double> Cf = {0.1,0.2,0.3,0.4,0.5};
  string method = "RungeKutta3";
  double cf = 0.5;

  // Setup sinewave data
  string sinefilename;
  vector<double> sinx_temp;
  sinx_temp.resize(extents[0]);

  sinefilename = "data_files/Sinewave-cf-"+to_string(cf)+".dat";
  ofstream sinefile;
  sinefile.open(sinefilename);
  for(int t=0;t<time_steps;t++)
  {
    if(t==0)
    {
      for(int i=0;i<sinx.size();i++)
      {
        sinefile << sinx[i] << " ";
      }
      sinefile << "\n";
    }
    for(int i=0;i<sinx.size();i++)
    {
      sinx_temp[i]=sin(i*(limits[0][1]-limits[0][0])/(static_cast<double>(extents[0]))-t*dx*cf);
      sinefile << sinx_temp[i] << " ";
    }
    sinefile << "\n";
  }
  sinefile.close();

  // Start main program
  U.set_all_data(sinx);
  U.update_ghostzone();

  // DO NOT USE thirdOrderDownstream!!!!
  // DOES NOT WORK

  TStepper(U,method,thirdOrderUpstream,P,time_steps,cf,cs,write_datafile);

}

