#include <iostream>
#include <vector>
using namespace std;

#include "domain.hpp"
#include "datamesh.cpp"

int main()
{
	vector<int> extents = {5,7,5};
	Mesh mesh = Mesh(extents.size(),extents);

	char test_mesh;
	char test_datamesh;
	char test_patch;

	cout << "True && True:" << (true && true) << endl;

	cout << "Test Mesh?: (T/F)" << endl;
	cin >> test_mesh;

	cout << "Test DataMesh?: (T/F)" << endl;
	cin >> test_datamesh;

	cout << "Test Patch?: (T/F)" << endl;
	cin >> test_patch;

	cout << "Extents are ";
       for(int i=0;i<extents.size();i++)
       {
		cout << extents[i] << " "; 
       }
       cout<<endl;
	
	if (test_mesh == 'T')
	{	
		cout << "Dimensions are " << mesh.get_dim() << endl;
		vector<int> exts = mesh.get_exts();
		cout << "Extents are ";
		for (int i=0; i<exts.size();i++){
			cout << exts[i] << ", ";
		}
		cout << endl;
		cout << "Total points are " << mesh.get_total_points() << endl;
	}

	if (test_datamesh == 'T')
	{
		DataMesh<int> dmA = DataMesh<int>(extents.size(),extents);
		DataMesh<int> dmB = DataMesh<int>(extents.size(),extents);
		DataMesh<bool> dmBool = DataMesh<bool>(extents.size(),extents);
		DataMesh<bool> dmBoolAlt = DataMesh<bool>(extents.size(),extents);
		cout << "DataMesh created" << endl;
		int coord;
		for(int i=0;i<extents[0];i++)
		{
			for(int j=0;j<extents[1];j++)
			{
				for(int k=0;k<extents[2];k++)
				{
					coord = k*(extents[0]*extents[1])+j*extents[0]+i;
					dmA.set_data_point(coord,i+j+k);
					dmB.set_data_point(coord,j);
					dmBool.set_data_point(coord,(bool) (k%2));
					dmBoolAlt.set_data_point(coord,(bool) (j%2));
				}
			}
		}
		cout << "Printing DataMesh A" << endl;
		dmA.print_datamesh();
		cout << "Printing DataMesh B" << endl;
		dmB.print_datamesh();
		dmA += dmB;
		cout << "Added B to A, and printing new A" << endl;
		dmA.print_datamesh();
		cout << "Making new DataMesh C by adding (new) A and B to it" << endl;
		DataMesh<int> dmC = dmA + dmB;
		dmC.print_datamesh();
		cout << "Doubling values in DM C by dmC*2" << endl;
		dmC *= 2;
		dmC.print_datamesh();

		cout << "Testing bool cases, Initial DataMesh dmBool:" << endl;
		dmBool.print_datamesh();
		cout << "Initial DataMesh dmBoolAlt:" << endl;
		dmBoolAlt.print_datamesh();
		dmBool += dmBoolAlt;
		cout << "+= dmBoolAlt to dmBool:" << endl;
		dmBool.print_datamesh();
	}
	

	if (test_patch == 'T')
	{
		cout << "Testing patch.cpp, domain of (-4,4),(-3,3),(-2,2)" << endl;
		cout << "Testing extents size: " << extents.size() << endl;

		vector<double> xlim = {-4.,4.};
		vector<double> ylim = {-3.,3.};
		vector<double> zlim = {-2.,2.};

		vector<vector<double>> limits = {xlim,ylim,zlim};

		Patch test_patch = Patch(extents.size(),extents,limits);
		test_patch.print();
	}




}
