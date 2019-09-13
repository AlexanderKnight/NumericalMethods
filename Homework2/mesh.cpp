#include <vector>
#include <iostream>
using namespace std;


#include "domain.hpp"

Mesh::Mesh(int dimensions, vector<int> extents)
/*
  Mesh Class Constructor
*/
{
	dim=dimensions;
	ext=extents;

	tot_points=1;
	for(int i=0;i<dim;i++)
	{
		tot_points *= ext[i];
	}
}

int Mesh::get_dim(void)
/*
  Gets dimension of Mesh
*/
{
	return dim;
}

vector<int> Mesh::get_exts(void)
/*
  Gets extents of Mesh
*/
{
	return ext;
}

int Mesh::get_total_points(void)
/*
  Gets total number of points in Mesh.
  Total points = Extents[0]*Extents[1]*....
*/
{
	return tot_points;
}
