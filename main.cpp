#include "SteinerTree.cpp"

int main(int argc, char **argv){
	srand(time(NULL));
	string file = string(argv[1]);
	SteinerTreeProblem STP(file);
	SteinerTree Individual;
	Individual.SteinerTreeproblem = &STP;
	Individual.restart();
}
