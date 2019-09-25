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
  stencil.push_back(tot_points);
	for(int i=0;i<dim;i++)
	{
		tot_points *= ext[i];
    if(i<(dim-1))
    {
      stencil.push_back(tot_points);
    }
	}
}

Mesh::Mesh(int dimensions, vector<int> extents, vector<int> gz_extents)
/*
  Mesh Class Constructor
*/
{
	dim=dimensions;
  gz_extents=gz_extents;
  sub_ext = extents;
  vector<int> total_extents;
  for(int i=0;i<dimensions;i++)
  {
    total_extents.push_back(extents[i]+2*gz_extents[i]);
  }
	ext=total_extents;

	tot_points=1;
  stencil.push_back(tot_points);
	for(int i=0;i<dim;i++)
	{
		tot_points *= total_extents[i];
    if(i<(dim-1))
    {
      stencil.push_back(tot_points);
    }
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

vector<int> Mesh::get_coords(int coord)
/*
 * Returns the integer location of a given coordinate
 */
{
	int remain_point;
  vector<int> coordinate;
  coordinate.resize(dim);
  for(int d=dim-1;d>=0;d--){
    if(d==dim-1)
    {
      coordinate[d] = (coord/stencil[d]);
      remain_point=coord%stencil[d];
    }
    else
    {
      coordinate[d]=(remain_point/stencil[d]);
      remain_point=remain_point%stencil[d];
    }
    //coordinates.set_data_point(i,coordinate);
  }
  return coordinate;
}
