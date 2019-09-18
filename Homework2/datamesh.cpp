#include <vector>
#include <iostream>
using namespace std;

#include "domain.hpp"


template<typename T>  
DataMesh<T>::DataMesh(int dimensions, vector<int> extents)
:Mesh(dimensions,extents)

/* Constructor for Class DataMesh
 */

{
	//Mesh(dimensions, extents);
	int product = 1;
	for(int i=0;i<extents.size();i++)
	{
		product *= extents[i];
	}
	data.resize(product);

}

template<typename T>  
int DataMesh<T>::get_total_points(void)

/* Gets total points in DataMesh.
 */

{
	return data.size();
}

template<typename T>  
void DataMesh<T>::set_all_data(vector<T> fill_data)

/* Given a vector of same size and data type, fills DataMesh with vector.
 * Useful to copy DataMesh.
 */

{
	if (fill_data.size() != data.size())
	{
		cout << "DataMesh and new data not same size" << endl;
	}
	else
	{
		for(int i=0;i<data.size();i++)
		{
			data[i] = fill_data[i];
		}
	}
}


template<typename T>  
vector<T> DataMesh<T>::get_all_data(void)

/* Gets all of data as one vector.
 */

{
	return data;
}

template<typename T>  
void DataMesh<T>::set_data_point(int coordinate, T data_point)
/* Sets individual point in DataMesh with given value.
 */
{
	data[coordinate]=data_point;
}

template<typename T>  
T DataMesh<T>::get_data_point(int coordinate)

/* Returns value from DataMesh coordinate.
 */

{
	return data[coordinate];
}

template<typename T>  
void DataMesh<T>::print_datamesh(void)

/* Prints the DataMesh to screen. Useful for debugging.
 *
 * WARNING: Currently assumes that DataMesh is 3D. 
 *          This will be fixed soon.
 */

{
	int coord;
	for(int i=0;i<ext[0];i++)
	{
		for(int j=0;j<ext[1];j++)
		{
			for(int k=0;k<ext[2];k++)
			{
				coord = k*(ext[0]*ext[1])+j*ext[0]+i;
				cout << data[coord] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
}

/*template<>  
void DataMesh<bool>::operator +=(DataMesh<bool> B)

/* Overloads the += operator for DataMesh<bool>.
 * Takes each point and does an AND operation for the pairs of values
 *   at the equivalent coordinates.
 *
 * Needs to be separated from the <int> and <double> += due to different
 *   handling for bools. (If we wanted an OR operation instead, then
 *   it would be the same)
*/

/*{
	assert(data.size() == B.get_total_points());
	assert(typeid(data) == typeid(B.get_all_data()));
	for(int i=0;i<data.size();i++)
	{
		data[i] = (data[i] && B.get_data_point(i));
	}
}*/

template<typename T>  
void DataMesh<T>::operator +=(DataMesh<T> B)

/* Overloads the += operator for DataMesh<int>, <double>, and <bool>.
 * Bool addition treated as AND operator.
 * Adds together each value at the equivalent coordinate.
 */

{
	assert(data.size() == B.get_total_points());
	assert(typeid(data) == typeid(B.get_all_data()));
  if(typeid(data[0])==typeid(bool))
  {
    for(int i=0;i<data.size();i++)
    {
      data[i] = (data[i] && B.get_data_point(i)); 
    }
  }
  else
  {
    for(int i=0;i<data.size();i++)
    {
      data[i] = (T)(data[i]+B.get_data_point(i));
    }
  }
}


/*template<>  
DataMesh<bool> DataMesh<bool>::operator +(DataMesh<bool> B)

/* Overloads the + operator to add two DataMesh<bool> together,
 * with the + operator being AND for the bool data at equivalent
 * coordinates.
 */

/*{
	assert(data.size() == B.get_total_points());
	assert(B.get_dim() == dim);
	assert(typeid(data) == typeid(B.get_all_data()));

	for(int i=0;i<B.get_dim();i++)
	{
		assert(B.get_exts()[i] == ext[i]);
	}
	
	DataMesh<bool> dm = DataMesh(B.get_dim(), B.get_exts());
	bool val;
	for(int i=0;i<B.get_total_points();i++)
	{
		val = (data[i] && B.get_data_point(i));
		dm.set_data_point(i,val);
	}
	return dm;
}*/

template<typename T>  
DataMesh<T> DataMesh<T>::operator +(DataMesh<T> B)

/* Overloads the + operator for <int> and <double> DataMesh's.
 * Adds data at equivalent coordinates
 */

{
	assert(data.size() == B.get_total_points());
  //assert(typeid(data[0]) == typeid(B.get_data_point(0)))
	assert(B.get_dim() == dim);

	for(int i=0;i<B.get_dim();i++)
	{
		assert(B.get_exts()[i] == ext[i]);
	}
	
	DataMesh<T> dm = DataMesh(B.get_dim(), B.get_exts());
	T val;
  if(typeid(data[0])==typeid(bool))
  {
    for(int i=0;i<B.get_total_points();i++)
    {
      val = (data[i] && B.get_data_point(i));
      dm.set_data_point(i,val);
    }
  }
  else
  {  
    for(int i=0;i<B.get_total_points();i++)
    {
      val = (T)(data[i]+B.get_data_point(i));
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
	for(int i=0;i<data.size();i++)
	{
		data[i] *= a;
	}
}
