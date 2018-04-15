main: MA.o SteinerTree.o utils.o main.cpp 
	g++ -o $@ main.cpp MA.o SteinerTree.o utils.o

SteinerTree.o: SteinerTree.h SteinerTree.cpp
	g++ -c -o $@ SteinerTree.cpp

MA.o: MA.h SteinerTree.h MA.cpp
	g++ -c -o $@ MA.cpp

utils.o: utils.cpp utils.h
	g++ -c -o $@ utils.cpp

clean:
	rm -f main SteinerTree.o MA.o utils.o
