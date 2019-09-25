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
		Mesh(int dim, vector<int> extents);
		Mesh(int dim, vector<int> extents, vector<int> gz_extents);
		int get_dim(void);
		vector<int> get_exts(void);
		int get_total_points(void);
    vector<int> get_coords(int coord);

	protected:
		int dim;  // dimensions
		vector<int> ext; // points along each axis
    vector<int> gz_extents;
    vector<int> sub_ext;
		int tot_points; // total points in mesh
    vector<int> stencil;
		
};

//------------------------------------------------------
//------------------------------------------------------

template <class T> class DataMesh : public Mesh
{
	public:
		DataMesh(int dimensions, vector<int> extents);
		DataMesh(int dimensions, vector<int> extents, vector<int> gz_extents, bool is_periodic=false);
    void setup_ghostzone(vector<int> extents,vector<int> gz_extents, bool is_periodic);
		void set_all_data(vector<T> fill_data);
    void update_gz(void);
		vector<T> get_all_data(void);
		void set_data_point(int coordinate, T data);
		T get_data_point(int coordinate);
		int get_total_points(void);
		void print(void);
    void print_ghostzone(void);
    void update_ghostzone(void);
    bool ghostzone(int coord);
		void operator +=(DataMesh<T> B);
		DataMesh<T> operator +(DataMesh<T> B);
		void operator *=(T a);
    void operator =(DataMesh<T> B);
    void write(ofstream& filename);
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
		Patch(int dimensions, vector<int>extents,vector<vector<double>> limits);
		void print(void);
		vector<double> get_coord(int coord);
    double dx(int);
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
    ComputeRHS(Patch patch, DataMesh<double> &U, DataMesh<double> &dU, double cs,
                string method, string differencing="centered");
    void ForwardEuler(Patch patch, DataMesh<double> &U, DataMesh<double> &dU, double cs,
                string differencing="centered");
};
#endif
