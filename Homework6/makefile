all:
	g++ -std=c++17 -o hw6 main.cpp domain.hpp mesh.cpp datamesh.cpp patch.cpp computerhs.cpp tstepper.cpp rootfinder.cpp rf.hpp

clean:
	rm hw6 data_files/*; find . -maxdepth 5 -type f -name "*png" -exec rm "{}" \;
