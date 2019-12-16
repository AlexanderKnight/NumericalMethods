all:
	g++ -std=c++17 -o hw5 main.cpp domain.hpp mesh.cpp datamesh.cpp patch.cpp computerhs.cpp tstepper.cpp
	g++ -std=c++17 -o rf rftest.cpp rf.hpp rootfinder.cpp 

rf:
	g++ -std=c++17 -o rf rftest.cpp rf.hpp rootfinder.cpp 

clean:
	rm hw5 rf data_files/*; find . -maxdepth 5 -type f -name "*png" -exec rm "{}" \;
