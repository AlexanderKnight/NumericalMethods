#ifndef DOMAIN_H
#define DOMAIN_H
#include <vector>
#include <assert.h>
#include <typeinfo>
#include <iostream>
#include <string>
using namespace std;


class Mesh
{
	public:
		Mesh(int dim, vector<int> &extents);
		Mesh(int dim, vector<int> &extents, vector<int> &gz_extents);
		int get_dim(void) const;
		vector<int> get_exts(void) const;
    vector<int> get_gz_exts(void) const;
		int get_total_points(void) const;
    vector<int> get_coords(int coord) const;

	protected:
		int dim;  // dimensions
		vector<int> ext; // points along each axis
    vector<int> gz_exts;
    vector<int> sub_ext;
		int tot_points; // total points in mesh
    vector<int> stencil;
		
};

//------------------------------------------------------
//------------------------------------------------------

template <class T> class DataMesh : public Mesh
{
	public:
		DataMesh(int dimensions, vector<int> &extents);
		DataMesh(int dimensions, vector<int> &extents, vector<int> &gz_extents, bool &is_periodic=false);
    void setup_ghostzone(vector<int> &extents,vector<int> &gz_extents, bool &is_periodic);
		void set_all_data(vector<T> &fill_data);
    void update_gz(void);
		vector<T> get_all_data(void) const;
		void set_data_point(int coordinate, T &data);
		T get_data_point(int coordinate) const;
		int get_total_points(void) const;
		void print(void) const;
    void print_ghostzone(void) const;
    void update_ghostzone(void);
    bool ghostzone(int coord) const;
		void operator +=(const DataMesh<T> &B);
		DataMesh<T> operator +(DataMesh<T> &B);
		void operator *=(T a);
    void operator =(DataMesh<T> &B);
    const T& operator [](int i) const;
    void write(ofstream &filename) const;
    void clean(void);
	private:
		vector<T> mesh_data;
    vector<bool> gz;
    int dimensions;

};


//----------------------------------------------------
//----------------------------------------------------

class Patch : public Mesh
{
	public:
		Patch(int dimensions, vector<int> &extents,vector<vector<double>> &limits);
		void print(void) const;
		vector<double> get_coord(int coord) const;
    double dx(int i) const;
	protected:
		vector<vector<double>> limits;
   	DataMesh<vector<double>> coordinates;
		vector<int> exts;
};

//----------------------------------------------------
//----------------------------------------------------

class ComputeRHS 
{
  public:
    //ComputeRHS(const DataMesh<double> &U, DataMesh<double> &dU, const double &cs, const double &dx,
                //string method, const double &dt, string differencing="upstream");
    void ForwardEuler(double (*fun)(const DataMesh<double> &, const int &, const double &),
                      const DataMesh<double> &U, DataMesh<double> &dU,
                      const double &cs, const double &dx, const double &dt);
    void RungeKutta3(double (*fun)(const DataMesh<double> &, const int&, const double&),
                     const DataMesh<double> &U, DataMesh<double> &dU,
                     const double &cs, const double &dx, const double &dt);
  private:
};

//----------------------------------------------------
//----------------------------------------------------

class TStepper
{
  public:
    TStepper(DataMesh<double> &U, string method, 
              double (*diff_func)(const DataMesh<double> &, const int &, const double &),
              double &dx, int &time_steps, double &cf, double &cs,
              bool write_datafile);
};

//=============================================================================
// 1st Order Derivatives
//=============================================================================

const double inline
centered(const DataMesh<double> &U, const int &i, const double &cs)
{
  return (U[i+1]-U[i-1])/2.;
}

double inline
downstream(const DataMesh<double> &U, const int &i, const double &cs)
{
  return (U[i+1]-U[i]);
}

double inline
upstream(const DataMesh<double> &U, const int &i, const double &cs)
{
  return (U[i]-U[i-1]);
}


//=============================================================================
// 3rd Order Derivatives
//=============================================================================

double inline
thirdOrderUpstream(const DataMesh<double> &U, const int &i, const double &cs)
{
  return (1./6.)*cs*(2.*U[i+1]+3.*U[i]-6.*U[i-1]+U[i-2]);
}

double inline
thirdOrderDownstream(const DataMesh<double> &U, const int &i, const double &cs)
{
  return (1./6.)*cs*(-2.*U[i-1]-3.*U[i]+6.*U[i+1]-U[i+2]);
}

//=============================================================================
// 4th Order Derivatives
//=============================================================================

double inline
fourthOrderCentered(const DataMesh<double> &U, const int &i, const double &cs)
{
  return (1./12.)*cs*(-1.*U[i+2]+8.*U[i+1]-8.*U[i-1]+U[i-2]);
}

#endif
