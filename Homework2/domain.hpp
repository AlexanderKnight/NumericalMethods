#ifndef DOMAIN_H
#define DOMAIN_H
#include <vector>
#include <assert.h>
#include <typeinfo>

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
		void operator *=(T a);
		DataMesh<T> operator +(DataMesh<T> B);
	private:
		vector<T> data;

};

template<typename T> inline DataMesh<T>::DataMesh(int dimensions, vector<int> extents):Mesh(dimensions,extents)
{
	Mesh(dimensions, extents);
	int product = 1;
	for(int i=0;i<extents.size();i++)
	{
		product *= extents[i];
	}
	data.resize(product);

}

template<typename T> inline void DataMesh<T>::set_all_data(vector<T> fill_data)
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

template<typename T> inline int DataMesh<T>::get_total_points(void)
{
	return data.size();
}

template<typename T> inline vector<T> DataMesh<T>::get_all_data(void)
{
	return data;
}

template<typename T> inline void DataMesh<T>::set_data_point(int coordinate, T data_point)
{
	data[coordinate]=data_point;
}

template<typename T> inline T DataMesh<T>::get_data_point(int coordinate)
{
	return data[coordinate];
}

template<typename T> inline void DataMesh<T>::print_datamesh(void)
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

template<> inline void DataMesh<bool>::operator +=(DataMesh<bool> B)
{
	assert(data.size() == B.get_total_points());
	assert(typeid(data) == typeid(B.get_all_data()));
	for(int i=0;i<data.size();i++)
	{
		data[i] = (data[i] && B.get_data_point(i));
	}
}

template<typename T> inline void DataMesh<T>::operator +=(DataMesh<T> B)
{
	assert(data.size() == B.get_total_points());
	assert(typeid(data) == typeid(B.get_all_data()));
	for(int i=0;i<data.size();i++)
	{
		data[i] += B.get_data_point(i);
	}
}


template<> inline DataMesh<bool> DataMesh<bool>::operator +(DataMesh<bool> B)
{
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
}

template<typename T> inline DataMesh<T> DataMesh<T>::operator +(DataMesh<T> B)
{
	assert(data.size() == B.get_total_points());
	assert(B.get_dim() == dim);

	for(int i=0;i<B.get_dim();i++)
	{
		assert(B.get_exts()[i] == ext[i]);
	}
	
	DataMesh<T> dm = DataMesh(B.get_dim(), B.get_exts());
	T val;
	for(int i=0;i<B.get_total_points();i++)
	{
		val = data[i]+B.get_data_point(i);
		dm.set_data_point(i,val);
	}
	return dm;
}

template<typename T> inline void DataMesh<T>::operator *= (T a)
{
	for(int i=0;i<data.size();i++)
	{
		data[i] *= a;
	}
}
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
	       	vector<vector<double>> coordinates;
		vector<int> exts;
};

#endif
