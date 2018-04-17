#ifndef __MA_H__
#define __MA_H__

#include "SteinerTree.h"

struct ExtendedIndividual {
	SteinerTree ind;
	int dist;
};

class MA {
	public:
		MA(int N_, double pc_, double pm_, double finalTime_);
		void run();
	private:
		//Parameters of MA
		int N;//Population Size
		double pc;//crossover probability
		double pm;//mutation probability
		double finalTime;//Seconds

		//Basic procedures of MA
		void initPopulation();
		void initDI();
		void selectParents();
		void crossover();
		void mutation();
		void localSearch();
		void replacement();

		//Internal attributes of MA
		vector< ExtendedIndividual * > population; 
		vector< ExtendedIndividual * > parents;
		vector< ExtendedIndividual * > offspring;
		double initialTime;
		double DI;
};

#endif
