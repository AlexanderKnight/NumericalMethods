#include <vector>
#include <iostream>
using namespace std;


#include "domain.hpp"

Mesh::Mesh(int dimensions, vector<int> extents)
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
