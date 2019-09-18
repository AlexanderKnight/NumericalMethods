output: main.o patch.o mesh.o
	g++ -std=c++17 -Wall -o hw2 main.o mesh.o patch.o

main.o: main.cpp domain.hpp
	g++ -c main.cpp
	
mesh.o: mesh.cpp domain.hpp
	g++ -c mesh.cpp

datamesh.o: datamesh.cpp domain.hpp
	g++ -c datamesh.cpp

patch.o: patch.cpp domain.hpp
	g++ -c patch.cpp

all:
	g++ -std=c++17 -Wall -o hw2 main.cpp mesh.cpp datamesh.cpp patch.cpp domain.hpp

clean:
	rm *.o hw2
