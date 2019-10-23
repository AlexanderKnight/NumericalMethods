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
  vector<vector<int> > extents = {{200}};
  //vector<vector<int> > extents = {{50}};
  vector<int> gz_extents = {2};
  vector<vector<double> > limits = {{0.,2.*M_PI}};
  bool is_periodic=true;
  double cs=1.;
  int time_steps;
  bool write_datafile = true;
  string method = "LLF";
  double cf = 0.1;
  vector<int> time_step = {4000};
  //vector<int> time_step = {1};

  for(int e=0;e<extents.size();e++)
  {
    cout << "Extents are " << extents[e][0] << endl;
    time_steps = time_step[e];
    DataMesh<double> U = DataMesh<double>(dimensions, extents[e],gz_extents,is_periodic);
    Patch P = Patch(dimensions, extents[e], limits);
    double dx = (2.*M_PI)/static_cast<double>(extents[e][0]);
    cout << "dx is " << dx << endl;
    
    // Sets up comparison gaussian wave
    vector<double> gaussx;
    gaussx.resize(extents[e][0]);
    double sigma = 0.2;
    double sigmaSq = pow(sigma,2);
    double coeff=pow(2.*M_PI*sigmaSq,-0.5);
    double num, exponential;
    double denom=2.*sigmaSq;
    cout << "Denom is " << denom << endl;
    
    for(int i=0;i<gaussx.size();i++)
    {
      num = i*(limits[0][1]-limits[0][0])/(static_cast<double>(extents[e][0]));
      exponential = exp(-1.*pow(num-M_PI,2)/(denom));
      cout << exponential << endl;
      gaussx[i]=coeff*exponential;
    }
    /*string sinefilename;
    vector<double> sinx_temp;
    sinx_temp.resize(extents[e][0]);
    sinefilename = "data_files/Sinewave-cf-"+to_string(cf)+"-"+to_string(extents[e][0])+".dat";
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
        sinx_temp[i]=sin(i*(limits[0][1]-limits[0][0])/(static_cast<double>(extents[e][0]))-(t+1)*dx*cf);
        sinefile << sinx_temp[i] << " ";
      }
      sinefile << "\n";
    }
    sinefile.close();
    */

    // Start main program
    U.set_all_data(gaussx);
    U.update_ghostzone();

    // Spacial derivative function is defined in header file
    // to ensure that all files have access to it.
    // DO NOT USE thirdOrderDownstream!!!!
    // DOES NOT WORK

    TStepper(U,method,thirdOrderUpstream,burger,dx,time_steps,cf,cs,write_datafile);
  }
}

