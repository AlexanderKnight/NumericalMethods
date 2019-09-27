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
		void operator +=(DataMesh<T> &B);
		DataMesh<T> operator +(DataMesh<T> &B);
		void operator *=(T &a);
    void operator =(DataMesh<T> &B);
    T operator [](int i) const;
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
    ComputeRHS(DataMesh<double> &U, DataMesh<double> &dU, double &cs, double &dx,
                string method, double &dt, string differencing="centered");
  private:
    vector<vector<double>> k;
    void ForwardEuler(DataMesh<double> &U, DataMesh<double> &dU, double &cs, double &dx, double &dt,
                string differencing="centered");
    void RungeKutta3(const DataMesh<double> &U, DataMesh<double> &dU, const double &cs,
                const double &dx, const double &dt);
};

//----------------------------------------------------
//----------------------------------------------------

class TStepper
{
  public:
    TStepper(DataMesh<double> &U, string method, string differencing,
              Patch &patch, int &time_steps, double &cf, double &cs,
              bool write_datafile);
};


#endif
