#include <vector>
#include <iostream>
using namespace std;

#include "domain.hpp"

Patch::Patch(int dimensions, vector<int>extents,
		vector<vector<double> > limits):Mesh(dimensions,extents)
/*
  Patch Class Constructor
*/

{
	//coordinates.resize(dimensions);

	vector<int> stencil = {1};

  //DataMesh<double> temp_dm = DataMesh<double>(dimensions,extents);
	for(int i=1;i<dimensions;i++)
	{
		stencil.push_back(stencil[i-1]*extents[i-1]);
		coordinates.push_back(DataMesh<double>(dimensions,extents));
	}
	
	int remain_point;
	for(int i=0;i<coordinates.size();i++)
	{
		for(int dim=dimensions-1;dim>=0;dim--){
			if(dim==dimensions-1)
			{
				coordinates[dim].set_data_point(i,limits[dim][0]
          +(i/stencil[dim])*(limits[dim][1]-limits[dim][0])/(extents[dim]-1));
				remain_point=i%stencil[dim];
			}
			else
			{
				coordinates[dim].set_data_point(i,limits[dim][0]
            +(remain_point/stencil[dim])
              *(limits[dim][1]-limits[dim][0])/(extents[dim]-1));
				remain_point=remain_point%stencil[dim];
			}
		}
	}
}

vector<double> Patch::get_coord(int coord)
/* 
  Given point, passes back spacial coordinates for that point
*/
{
	vector<double> pass_coord;
	for(int i=0;i<coordinates.size();i++)
	{
		pass_coord.push_back(coordinates[i].get_data_point(coord));
	}
	return pass_coord;
}

void Patch::print(void)
/*
 *  Prints the patch coordinates in a grid pattern, for easier debugging.
 *  Currently requires cases for 1D, 2D, and 3D. Looking to improve to a
 *  general case.
 *  
 *  Format Example: 2x2x2 patch (x,y,z)
 *
 *  (0,0,0) (1,0,0)
 *  (0,1,0) (1,1,0)
 *
 *  (0,0,1) (1,0,1)
 *  (0,1,1) (1,1,1)
 *
 *  Column: costant x,
 *  Row: constant y,
 *  Block: constant z
 */
{

	int coor;
	if(ext.size()==3)
	{
		for(int k=0;k<this->ext[2];k++)
		{
			for(int j=0;j<this->ext[1];j++)
			{
				for(int i=0;i<this->ext[0];i++)
				{
					coor = k*this->ext[0]*this->ext[1]+j*this->ext[0]+i;
					cout << "(" << coordinates[0].get_data_point(coor) << "," 
						<< coordinates[1].get_data_point(coor) << "," 
						<< coordinates[2].get_data_point(coor) << ") ";
				}
				cout << endl;
			}
			cout << endl;
		}
	}
	else if (this->ext.size()==2)
	{
		for(int j=0;j<this->ext[1];j++)
		{
			for(int i=0;i<this->ext[0];i++)
			{
				coor = j*this->ext[0]+i;
				cout<<"("<<coordinates[0].get_data_point(coor)<<","
              <<coordinates[1].get_data_point(coor)<<") ";
			}
			cout << endl;
		}
	}
	else
	{
		for(int i=0;i<this->tot_points;i++)
		{
			cout << coordinates[0].get_data_point(i) << " ";
		}
	}
}



