#include <vector>
#include <iostream>
#include <typeinfo>
using namespace std;

#include "domain.hpp"
#include "datamesh.cpp"

Patch::Patch(int dimensions, vector<int>extents,
		vector<vector<double> > limits):Mesh(dimensions,extents),
  coordinates(DataMesh<vector<double>>(dimensions,extents))
/*
  Patch Class Constructor
*/

{
// Reminder that coordinates is DataMesh<vector<double>>
	vector<int> stencil = {1};
  //coordinates=DataMesh<vector<double>>(dimensions,extents);
	for(int i=1;i<dimensions;i++)
	{
		stencil.push_back(stencil[i-1]*extents[i-1]);
cout << "Stencil is " << stencil[i-1]*extents[i-1] << endl;
	}
	
//----------------------------------------------
	int remain_point;
  vector<double> val= {0.,0.,0.};
	for(int i=0;i<coordinates.get_total_points();i++)
	{
		for(int dim=dimensions-1;dim>=0;dim--){
			if(dim==dimensions-1)
			{
        val[dim] = limits[dim][0]+(i/stencil[dim])*(limits[dim][1]-limits[dim][0])/(extents[dim]-1);
				//coordinates[dim].set_data_point(i, val);
				remain_point=i%stencil[dim];
			}
			else
			{
        val[dim]=limits[dim][0]+(remain_point/stencil[dim])
              *(limits[dim][1]-limits[dim][0])/(extents[dim]-1);
				//coordinates[dim].set_data_point(i,val);
				remain_point=remain_point%stencil[dim];
			}
      coordinates.set_data_point(i,val);
		}
	}
}

vector<double> Patch::get_coord(int coord)
/* 
  Given point, passes back spacial coordinates for that point
*/
{
	return coordinates.get_data_point(coord);
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
					cout << "(" << coordinates.get_data_point(coor)[0] << "," 
						<< coordinates.get_data_point(coor)[1] << "," 
						<< coordinates.get_data_point(coor)[2] << ") ";
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
				cout<<"("<<coordinates.get_data_point(coor)[0]<<","
              <<coordinates.get_data_point(coor)[1]<<") ";
			}
			cout << endl;
		}
	}
	else
	{
		for(int i=0;i<this->tot_points;i++)
		{
			cout << coordinates.get_data_point(i)[0] << " ";
		}
	}
  cout << endl;
}

double Patch::dx(int dim)
{
  return coordinates.get_data_point(stencil[dim])[dim]
            -coordinates.get_data_point(0)[dim];
}
  



