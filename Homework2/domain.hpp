#ifndef DOMAIN_H
#define DOMAIN_H
#include <vector>
#include <assert.h>
#include <typeinfo>
using namespace std;

class Mesh
{
	public:
		Mesh(int dim, vector<int> extents);
		int get_dim(void);
		vector<int> get_exts(void);
		int get_total_points(void);

	protected:
		int dim;  // dimensions
		vector<int> ext; // points along each axis
		int tot_points; // total points in mesh
		
};

//------------------------------------------------------
//------------------------------------------------------

/*
    NOTE: Currently, all member functions of class 
        DataMesh are in this header. Long term planning is to 
        shift them to their own file, and just have the declarations here.
*/

template <class T> class DataMesh : public Mesh
{
	public:
		DataMesh(int dimensions, vector<int> extents);
		void set_all_data(vector<T> fill_data);
		vector<T> get_all_data(void);
		void set_data_point(int coordinate, T data);
		T get_data_point(int coordinate);
		int get_total_points(void);
		void print_datamesh(void);
		void operator +=(DataMesh<T> B);
		DataMesh<T> operator +(DataMesh<T> B);
		void operator *=(T a);
		//void operator +=(DataMesh<bool> B);
		//DataMesh<bool> operator +(DataMesh<bool> B);
	private:
		vector<T> data;

};


//----------------------------------------------------
//----------------------------------------------------

class Patch : public Mesh
{
	public:
		Patch(int dimensions, vector<int>extents,vector<vector<double>> limits);
		void print(void);
		vector<double> get_coord(int coord);
	protected:
		vector<vector<double>> limits;
   	vector<DataMesh<double>> coordinates;
		vector<int> exts;
};

#endif
