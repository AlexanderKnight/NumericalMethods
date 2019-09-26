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
  vector<int> extents = {300};
  vector<int> gz_extents = {1};
  vector<vector<double> > limits = {{0.,2.*M_PI}};
  bool is_periodic=true;
  double cs=1.;
  double cf;
  int time_steps = 300;
  
  vector<double> sinx;
  sinx.resize(extents[0]);
  for(int i=0;i<sinx.size();i++)
  {
    sinx[i]=sin(i*(limits[0][1]-limits[0][0])/(static_cast<double>(extents[0])));
  }
  DataMesh<double> U = DataMesh<double>(dimensions, extents,gz_extents,is_periodic);
  Patch P = Patch(dimensions, extents, limits);
  double dx = P.dx(0);
  vector<string> differencing = {"centered", "downstream", "upstream"};
  vector<double> Dt = {0.1,0.2,0.3,0.4,0.5};
  string filename;
  //double dt = 0.25*dx;
  //double dt = 0.00628947;
  cout << "dx is " << dx << endl;
  cout << "Total points are " << U.get_total_points() << endl;

  for(int diffs=0;diffs<differencing.size();diffs++)
  {
    cout << differencing[diffs] << endl;
    for(int dt=0;dt<Dt.size();dt++)
    {
      U.set_all_data(sinx);
      U.update_ghostzone();
      U.print();
      cout << "    " << Dt[dt] << endl;
      filename = "data_files/ForwardEuler-"+differencing[diffs]+
                      "-SineWave-dt-"+to_string(Dt[dt])+".dat";
      ofstream datafile;
      datafile.open(filename);
      DataMesh<double> dU = DataMesh<double>(dimensions,extents,gz_extents,is_periodic);
      U.write(datafile);
      for(int t=0;t<time_steps;t++)
      {
        if(t%50==0)
        {
          cout << "        " << t << endl;
        }
        ComputeRHS(P,U,dU,cs,"ForwardEuler",Dt[dt],differencing[diffs]);
        cf = dx*Dt[dt];
        dU *= cf;
        dU += U;
        U = dU;
        U.update_ghostzone();
        U.write(datafile);
      }
      datafile.close();
      U.clean();
    }
  }

  string sinefilename;
  vector<double> sinx_temp;
  sinx_temp.resize(extents[0]);
  for(int dt=0;dt<Dt.size();dt++)
  {
    sinefilename = "data_files/Sinewave-dt-"+to_string(Dt[dt])+".dat";
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
        sinx_temp[i]=sin(i*(limits[0][1]-limits[0][0])/(static_cast<double>(extents[0]))-t*dx*Dt[dt]);
        sinefile << sinx_temp[i] << " ";
      }
      sinefile << "\n";
    }
    sinefile.close();
  }
}

