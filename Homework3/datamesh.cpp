#include <vector>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "domain.hpp"


template<typename T>  
DataMesh<T>::DataMesh(int dimensions, vector<int> &extents)
:Mesh(dimensions,extents)

/* Constructor for Class DataMesh
 */

{
  dimensions = dimensions;
	int product = 1;
	for(int i=0;i<extents.size();i++)
	{
		product *= extents[i];
	}
	mesh_data.resize(product);

}

template<typename T>  
DataMesh<T>::DataMesh(int dimensions, vector<int> &extents, 
                      vector<int> &gz_extents, bool &is_periodic)
:Mesh(dimensions,extents,gz_extents)

/* Constructor for Class DataMesh with ghost zone
 */

{
  dimensions = dimensions;
	int product = 1;
  assert (extents.size() == gz_extents.size());
	for(int i=0;i<extents.size();i++)
	{
		product *= extents[i]+2*gz_extents[i];
	}
	mesh_data.resize(product);
  gz.resize(product);
  setup_ghostzone(ext,gz_extents,is_periodic);

}

template<typename T>
void DataMesh<T>::setup_ghostzone(vector<int> &extents,vector<int> &gz_extents, bool &is_periodic)
    
/* 
 * Builds ghost zone mask, covers edges of boundary of simulation
 */
{
  

  vector<int> coordinates;
  for(int coord=0;coord<mesh_data.size();coord++)
  {
    coordinates = get_coords(coord);
    for(int dim=0;dim<extents.size();dim++)
    {
      if((coordinates[dim] < gz_extents[dim]) ||
           (coordinates[dim] >= extents[dim]-gz_extents[dim]))
      {
        gz[coord]=1;
      }
    }
  }
  update_ghostzone();
}
     
template<typename T>
void DataMesh<T>::update_ghostzone(void)

/*
 * For a periodic boundary, this grabs points on the opposite side
 * of the domain to use for derivatives at the boundary
 */
{
  for(int dim=0;dim<ext.size();dim++)
  {
    for(int i=0;i<mesh_data.size();i++)
    {
      if(gz[i]==1)
      {
        if(get_coords(i)[dim]<ext[dim]/2)
        {
          if(gz[i+stencil[dim]*sub_ext[dim]]==0)
          {
            mesh_data[i]=mesh_data[i+stencil[dim]*sub_ext[dim]];
          }
        }
        else
        {
          if(gz[i-stencil[dim]*sub_ext[dim]]==0)
          {
            mesh_data[i]=mesh_data[i-stencil[dim]*sub_ext[dim]];
          }
        }
      }
    }
  }
}

template<typename T>
bool DataMesh<T>::ghostzone(int coord) const

/* Returns if coordinate is a ghost zone point or not
 */
{
  return gz[coord];
}

template<typename T>  
int DataMesh<T>::get_total_points(void) const

/* Gets total points in DataMesh.
 */

{
	return mesh_data.size();
}

template<typename T>  
void DataMesh<T>::set_all_data(vector<T> &fill_data)

/* Given a vector of same size and data type, fills DataMesh with vector.
 * Useful to copy DataMesh.
 */

{
  int k =0;
  if(gz.size()==0) //If there is no ghost zone, then can fill right away
  {
    if (fill_data.size() != mesh_data.size())
    {
      cout << "DataMesh and new data not same size" << endl;
    }
    else
    {
      for(int i=0;i<mesh_data.size();i++)
      {
        mesh_data[i] = fill_data[i];
      }
    }
  }
  /* If there is a ghost zone, then filling is tricky. 
   * Example: Let the DM be a NxM, with a 3 cell ghost zone.
   *          Then the actual size is (N+2*3)x(M+2*3). 
   *          If the user wants to fill it, they will use a NxM,
   *          the same size as they specified when it was created.
   *          This part skips ghost zone points and only changes the
   *          non-ghost zone points (of which is there is N*M of them).
   */
  else
  {
    for(int i=0;i<mesh_data.size();i++)
    {
      if(gz[i]==0)
      {
        mesh_data[i] = fill_data[i-k];
      }
      else
      {
        k++;
      }
    }
	}
}


template<typename T>  
vector<T> DataMesh<T>::get_all_data(void) const

/* Gets all of data as one vector.
 */

{
	return mesh_data;
}

template<typename T>  
void DataMesh<T>::set_data_point(int coordinate, T &data_point)

/* Sets individual point in DataMesh with given value.
 */
{
	mesh_data[coordinate]=data_point;
}

template<typename T>  
T DataMesh<T>::get_data_point(int coordinate) const

/* Returns value from DataMesh coordinate.
 */

{
	return mesh_data[coordinate];
}

template<typename T>  
void DataMesh<T>::print(void) const

/* Prints the DataMesh to screen. Useful for debugging.
 */

{
	int coor;
	if(ext.size()==3)
	{
		for(int k=0;k<ext[2];k++)
		{
			for(int j=0;j<ext[1];j++)
			{
				for(int i=0;i<ext[0];i++)
				{
          cout << mesh_data[k*ext[0]*ext[1]+j*ext[0]+i] << " ";
				}
				cout << endl;
			}
			cout << endl;
		}
	}
	else if (ext.size()==2)
	{
		for(int j=0;j<ext[1];j++)
		{
			for(int i=0;i<ext[0];i++)
			{
          cout << mesh_data[j*ext[0]+i] << " ";
			}
			cout << endl;
		}
	}
	else
	{
		for(int i=0;i<tot_points;i++)
		{
			cout << mesh_data[i] << " ";
		}
	}
  cout << endl;
}

template<typename T>  
void DataMesh<T>::print_ghostzone(void) const

/* Prints the DataMesh to screen. Useful for debugging.
 */

{
	if(ext.size()==3)
	{
		for(int k=0;k<ext[2];k++)
		{
			for(int j=0;j<ext[1];j++)
			{
				for(int i=0;i<ext[0];i++)
				{
          cout << gz[k*ext[0]*ext[1]+j*ext[0]+i] << " ";
				}
				cout << endl;
			}
			cout << endl;
		}
	}
	else if (ext.size()==2)
	{
		for(int j=0;j<ext[1];j++)
		{
			for(int i=0;i<ext[0];i++)
			{
          cout << gz[j*ext[0]+i] << " ";
			}
			cout << endl;
		}
	}
	else
	{
		for(int i=0;i<tot_points;i++)
		{
			cout << gz[i] << " ";
		}
	}
  cout << endl;
}

template<typename T>
void DataMesh<T>::write(ofstream &filename) const

/* Saves DM to file, specified by a stream. Currently saves entire DM in
 * one row. Future plans is to make that more workable for most plotting 
 * programs. As it is, numpy will import the data and plot correctly.
 *
 * Future: auto remove ghost zone points, so they don't have to be removed
 *          post processing
 */
{
	if(ext.size()==3)
	{
		for(int k=0;k<ext[2];k++)
		{
			for(int j=0;j<ext[1];j++)
			{
				for(int i=0;i<ext[0];i++)
				{
          filename << mesh_data[k*ext[0]*ext[1]+j*ext[0]+i] << " ";
				}
				filename << "\n";
			}
			filename << "\n";
		}
	}
	else if (ext.size()==2)
	{
		for(int j=0;j<ext[1];j++)
		{
			for(int i=0;i<ext[0];i++)
			{
          filename << mesh_data[j*ext[0]+i] << " ";
			}
			filename << "\n";
		}
	}
	else
	{
		for(int i=0;i<tot_points;i++)
		{
			filename << mesh_data[i] << " ";
		}
	}
  filename << "\n";
}
  

template<typename T>  
inline void DataMesh<T>::operator +=(const DataMesh<T> &B)

/* Overloads the += operator for DataMesh<int>, <double>, and <bool>.
 * Bool addition treated as AND operator.
 * Adds together each value at the equivalent coordinate.
 */

{
	assert(mesh_data.size() == B.get_total_points());
	assert(typeid(mesh_data[0]) == typeid(B.get_data_point(0)));
  if(typeid(T)==typeid(bool))
  {
    for(int i=0;i<mesh_data.size();i++)
    {
      mesh_data[i] = mesh_data[i]*B.get_data_point(i); 
    }
  }
  else
  {
    for(int i=0;i<mesh_data.size();i++)
    {
      mesh_data[i] = (T)(mesh_data[i]+B.get_data_point(i));
    }
  }
}

template<typename T>  
DataMesh<T> DataMesh<T>::operator +(DataMesh<T> &B)

/* Overloads the + operator for <int> and <double> DataMesh's.
 * Adds data at equivalent coordinates
 */

{
	assert(mesh_data.size() == B.get_total_points());
  //assert(typeid(data[0]) == typeid(B.get_data_point(0)))
	assert(B.get_dim() == dim);

	for(int i=0;i<B.get_dim();i++)
	{
		assert(B.get_exts()[i] == ext[i]);
	}
	
	DataMesh<T> dm = DataMesh(B.get_dim(), B.get_exts());
	T val;
  if(typeid(T)==typeid(bool))
  {
    for(int i=0;i<B.get_total_points();i++)
    {
      val = mesh_data[i]*B.get_data_point(i);
      dm.set_data_point(i,val);
    }
  }
  else
  {  
    for(int i=0;i<B.get_total_points();i++)
    {
      val = (T)(mesh_data[i]+B.get_data_point(i));
      dm.set_data_point(i,val);
    }
  }
	return dm;
}

template<typename T>  
void DataMesh<T>::operator *= (T a)

/* Overloads *= operator, for multiplying a <int> or <double> 
 * Datamesh by a scalar of the same type.
 */

{
  assert(typeid(a) != typeid(bool));
	for(int i=0;i<mesh_data.size();i++)
	{
		mesh_data[i] *= a;
	}
}

template<typename T>  
void DataMesh<T>::operator =(DataMesh<T> &B)

/* Overloads the + operator for <int> and <double> DataMesh's.
 * Adds data at equivalent coordinates
 */

{
	assert(mesh_data.size() == B.get_total_points());
	assert(B.get_dim() == dim);

	for(int i=0;i<B.get_dim();i++)
	{
		assert(B.get_exts()[i] == ext[i]);
	}
	
  for(int i=0;i<B.get_total_points();i++)
  {
    mesh_data[i]= B.get_data_point(i);
  }
}

template<typename T>
const T& DataMesh<T>::operator [](int i) const
{
  assert(i < mesh_data.size());
  return mesh_data[i];
}

template<typename T>
void DataMesh<T>::clean(void)

/* Resets all values in the DM to zero. 
 * For <int,double,bool> this is (0,0.,false).
 *
 * Created as sometimes the ghost zones can cause problems, so this just 
 * wipes everything.
 */
{
  if(typeid(T)==typeid(bool) || typeid(T)==typeid(int))
  {
    for(int i=0;i<mesh_data.size();i++)
    {
      mesh_data[i]=0;
    }
  }
  else
  {
    for(int i=0;i<mesh_data.size();i++)
    {
      mesh_data[i]=0.;
    }
  }
}
