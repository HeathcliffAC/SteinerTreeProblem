#include "MA.h"

int main(int argc, char **argv){
	int N = 50;
	double pc = 0.9;
	double pm = 0.01;
	double finalTime = 25 * 60;
	MA ma(N, pc, pm, finalTime);
	srand(time(NULL));
	//string file = string(argv[1]);
	SteinerTreeProblem STP;
	SteinerTree::SteinerTreeproblem = &STP;
	ma.run();
}
