#include <iostream>
#include <vector>
using namespace std;

#include "domain.hpp"

Patch::Patch(int dimensions, vector<int>extents,vector<vector<double> > limits):Mesh(dimensions,extents)
{
	coordinates.resize(this->tot_points);

	vector<int> stencil = {1};

	for(int i=1;i<dimensions;i++)
	{
		stencil.push_back(stencil[i-1]*extents[i-1]);
	}
	
	int remain_point;
	for(int i=0;i<coordinates.size();i++)
	{
		coordinates[i].resize(dimensions);
		for(int dim=dimensions-1;dim>=0;dim--){
			if(dim==dimensions-1)
			{
				//cout << "coord: " << i << endl;
				//cout << "dim: " << dim << endl;
				//cout << "stencil[dim]: " << stencil[dim] << endl;
				//cout << "coord/stencil[dim]: " << i/stencil[dim] << endl;
				//cout << limits[dim][1] << endl;
				//cout << limits[dim][0] << endl;
				//cout << extents[dim] << endl;
				//cout << limits[dim][0]/(extents[dim]-1) << endl;
				coordinates[i][dim]=limits[dim][0]+(i/stencil[dim])*(limits[dim][1]-limits[dim][0])/(extents[dim]-1);
				remain_point=i%stencil[dim];
			}
			else
			{
				coordinates[i][dim]=limits[dim][0]+(remain_point/stencil[dim])*(limits[dim][1]-limits[dim][0])/(extents[dim]-1);
				remain_point=remain_point%stencil[dim];
			}
		}
	}
}

vector<double> Patch::get_coord(int coord)
{
	return coordinates[coord];
}

void Patch::print(void)
{

	int coor;
	cout << "Extent size is " << this->ext.size() << endl;
	if(ext.size()==3)
	{
		for(int k=0;k<this->ext[2];k++)
		{
			for(int j=0;j<this->ext[1];j++)
			{
				for(int i=0;i<this->ext[0];i++)
				{
					coor = k*this->ext[0]*this->ext[1]+j*this->ext[0]+i;
					cout << "(" << coordinates[coor][0] << "," << coordinates[coor][1] << "," << coordinates[coor][2] << ") ";
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
				cout << "(" << coordinates[coor][0] << "," << coordinates[coor][1] << ") ";
			}
			cout << endl;
		}
	}
	else
	{
		for(int i=0;i<this->tot_points;i++)
		{
			cout << coordinates[i][0] << " ";
		}
	}
}



