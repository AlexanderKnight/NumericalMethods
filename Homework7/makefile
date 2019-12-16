all:
	g++ -std=c++17 -o hw7 main.cpp domain.hpp mesh.cpp datamesh.cpp patch.cpp computerhs.cpp tstepper.cpp rootfinder.cpp minimizer.cpp min.hpp rf.hpp

clean:
	rm hw* data_files/*; find . -maxdepth 5 -type f -name "*png" -exec rm "{}" \;
