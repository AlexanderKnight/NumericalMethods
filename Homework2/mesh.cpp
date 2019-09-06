#include <vector>
#include <list>
#include <iostream>
useing namespace std;

class Mesh
{
	private:
		int dim;  // dimensions
		vector<int> ext; // points along each axis
		int tot_points; // total points in mesh
		
	public:
		Mesh(int dim, vector<int> extents);
};

Mesh::Mesh(int dimensions, vector<int> extents)
{
	dim=dimensions;
	ext=extents;

	tot_points=1;
	for(int i=0;i<n;i++)
	{
		tot_points *= ext[i];
	}
}

int Mesh::get_dim(void)
{
	return dim;
}

vector<int> Mesh::get_exts(void)
{
	return ext;
}

int Mesh::get_total_points(void)
{
	return tot_points;
}
