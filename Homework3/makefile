all:
	g++ -std=c++17 -o hw3 main.cpp domain.hpp mesh.cpp datamesh.cpp patch.cpp computerhs.cpp tstepper.cpp

clean:
	rm hw3 data_files/*; find . -maxdepth 5 -type f -name "*png" -exec rm "{}" \;
