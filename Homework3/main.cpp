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
  // Setup simulation parameters
  int dimensions = 1;
  vector<vector<int> > extents = {{100},{200},{400},{800}};
  //vector<vector<int> > extents = {{50}};
  vector<int> gz_extents = {2};
  vector<vector<double> > limits = {{0.,2.*M_PI}};
  bool is_periodic=true;
  double cs=1.;
  int time_steps;
  bool write_datafile = true;
  string method = "RungeKutta3";
  double cf = 0.1;
  vector<int> time_step = {1000,2000,4000,8000};
  //vector<int> time_step = {1};

  for(int e=0;e<extents.size();e++)
  {
    cout << "Extents are " << extents[e][0] << endl;
    time_steps = time_step[e];
    //time_steps = (int)(extents[e][0]/cf);
    DataMesh<double> U = DataMesh<double>(dimensions, extents[e],gz_extents,is_periodic);
    Patch P = Patch(dimensions, extents[e], limits);
    //double dx = P.dx(0);
    double dx = (2.*M_PI)/static_cast<double>(extents[e][0]);
    cout << "dx is " << dx << endl;
    
    // Sets up comparison sine wave
    cout << "Setting up sine wave" << endl;
    vector<double> sinx;
    sinx.resize(extents[e][0]);
    for(int i=0;i<sinx.size();i++)
    {
      sinx[i]=sin(i*(limits[0][1]-limits[0][0])/(static_cast<double>(extents[e][0])));
    }
    string sinefilename;
    vector<double> sinx_temp;
    sinx_temp.resize(extents[e][0]);
    sinefilename = "data_files/Sinewave-cf-"+to_string(cf)+"-"+to_string(extents[e][0])+".dat";
    ofstream sinefile;
    sinefile.open(sinefilename);
    cout << "Updating Sine Wave" << endl;
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
        sinx_temp[i]=sin(i*(limits[0][1]-limits[0][0])/(static_cast<double>(extents[e][0]))-(t+1)*dx*cf);
        sinefile << sinx_temp[i] << " ";
      }
      sinefile << "\n";
    }
    sinefile.close();

    // Start main program
    cout << "Setting sine data to U" << endl;
    U.set_all_data(sinx);
    U.update_ghostzone();

    // DO NOT USE thirdOrderDownstream!!!!
    // DOES NOT WORK

    cout << "Starting TStepper" << endl;
    TStepper(U,method,thirdOrderUpstream,dx,time_steps,cf,cs,write_datafile);
    cout << "TStepper done" << endl;
  }
}

