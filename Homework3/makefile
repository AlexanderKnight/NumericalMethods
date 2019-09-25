#output: test.o patch.o mesh.o
#	g++ -std=c++17 -Wall -o hw2 test.o mesh.o patch.o

#test.o: test.cpp domain.hpp
#	g++ -c test.cpp
	
#mesh.o: mesh.cpp domain.hpp
#	g++ -c mesh.cpp

#datamesh.o: datamesh.cpp domain.hpp
#	g++ -c datamesh.cpp

#patch.o: patch.cpp domain.hpp
#	g++ -c patch.cpp

all:
	g++ -std=c++17 -o hw3 main.cpp domain.hpp mesh.cpp datamesh.cpp patch.cpp computerhs.cpp

clean:
	rm hw3 *.dat data_files/*; find . -maxdepth 5 -type f -name "*png" -exec rm "{}" \;
