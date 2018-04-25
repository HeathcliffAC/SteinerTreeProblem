#include <sys/time.h>
#include <iostream>
#include <signal.h>

#include "MA.h"
#include "utils.h"

using namespace std;

void printer(int signal){
	finished = true;
}

MA::MA(int N_, double pc_, double pm_, double finalTime_){
	signal(SIGTERM, printer);
	//signal(SIGALRM, printer);
	N = N_;
	if (N % 2){ cerr << "El tam. de poblacion debe ser par" << endl; exit(-1); }
	pc = pc_;
	pm = pm_;
	finalTime = finalTime_;
	struct timeval currentTime; 
	gettimeofday(&currentTime, NULL);
	initialTime = (double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6;
}

void MA::initPopulation(){
	struct timeval initTime; 
	gettimeofday(&initTime, NULL);
	double time = ((double) (initTime.tv_sec) * 1.0e6 + (double) (initTime.tv_usec))/1.0e6;
	for (int i = 0; i < N; i++){
		//cout << "Crea ind " << i << endl;
		ExtendedIndividual *ei = new ExtendedIndividual();
		ei->ind.restart();
		population.push_back(ei);
		if ((population.size() >= 4) && (population.size() % 2 == 0)){
			struct timeval currentTime;
			gettimeofday(&currentTime, NULL);
			double time2 = ((double) (currentTime.tv_sec) * 1.0e6 + (double) (currentTime.tv_usec))/1.0e6;
			double elapsed = time2 - time;
			//cout << "Va: " << elapsed << endl;
			if (elapsed > 30 * 60.0 / 2500){
				N = i + 1;
			}
		}
	}
	//cout << "Tam. pob: " << N << endl;
}

//Select parents with binary selection
void MA::selectParents(){
	parents.clear();
	if (generation > 100){//Mutation is active
		parents.push_back(population[0]);
		parents.push_back(population[0]);
	}
	for (int i = parents.size(); i < N; i++){
		int first = getRandomInteger0_N(N - 1);
		int second = getRandomInteger0_N(N - 1);
		if (population[first]->ind.fitness <= population[second]->ind.fitness){
			parents.push_back(population[first]);
		} else {
			parents.push_back(population[second]);
		}
	}
}

void MA::crossover(){
	for (int i = 0; i < parents.size(); i++){
		ExtendedIndividual *ei = new ExtendedIndividual();
		*ei = *parents[i];
		offspring.push_back(ei);
	}
	for (int i = 0; i < offspring.size(); i+=2){
		if (generateRandomDouble0_Max(1) <= pc){
			offspring[i]->ind.dependentCrossover(offspring[i+1]->ind);
		}
	}
}

void MA::mutation(){
	if (generation > 100){
		for (int i = 0; i < offspring.size(); i++){
			offspring[i]->ind.dependentMutation(pm);
		}
	}
}

void MA::localSearch(){
	for (int i = 0; i < offspring.size(); i++){
		offspring[i]->ind.localSearch();
	}
}


void MA::replacement(){
	vector < ExtendedIndividual* > all;
	
	//Join population and offspring
	for (int i = 0; i < population.size(); i++){
		//cout << "Valor: " << population[i]->ind.fitness << endl;
		all.push_back(population[i]);
		all.back()->dist = INT_MAX;
	}
	population.clear();
	for (int i = 0; i < offspring.size(); i++){
		//cout << "Valor: " << offspring[i]->ind.fitness << endl;
		all.push_back(offspring[i]);
		all.back()->dist = INT_MAX;
	}
	offspring.clear();
	
	//Select best solution
	int indexBest = 0;
	for (int i = 1; i < all.size(); i++){
		if (all[i]->ind.fitness < all[indexBest]->ind.fitness){
			indexBest = i;
		}
	}
	//cout << "Mete a: " << all[indexBest]->ind.fitness << endl;
	population.push_back(all[indexBest]);
	all[indexBest] = all.back();
	all.pop_back();

	struct timeval currentTime; 
	gettimeofday(&currentTime, NULL);
	double elapsedTime = (double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6;
	elapsedTime -= initialTime;

	//Select next N - 1 solution
	double D = DI - DI * elapsedTime / finalTime;
	//cout << "Distancia requerida: " << D << endl;
	set<long long> acceptedFitness;
	acceptedFitness.insert(population[0]->ind.fitness);
	while(population.size() != N){
		//Update distances
		for (int i = 0; i < all.size(); i++){
			all[i]->dist = min(all[i]->dist, all[i]->ind.getDistance(population.back()->ind));
		}
		//Select best option
		indexBest = 0;
		for (int i = 1; i < all.size(); i++){
			bool betterInDist =	(all[i]->dist > all[indexBest]->dist);
			bool eqInDist = (all[i]->dist == all[indexBest]->dist);
			bool betterInFit = (all[i]->ind.fitness < all[indexBest]->ind.fitness);
			bool eqInFit = (all[i]->ind.fitness == all[indexBest]->ind.fitness);
			/*if ((acceptedFitness.count(all[indexBest]->ind.fitness) != 0) && (acceptedFitness.count(all[i]->ind.fitness) == 0)){
				indexBest = i;
			} else if ((acceptedFitness.count(all[indexBest]->ind.fitness) == 0) && (acceptedFitness.count(all[i]->ind.fitness) != 0)){
			} else */
			if (all[indexBest]->dist < D){//Do not fulfill distance requirement
				if ((betterInDist) || (eqInDist && betterInFit)){
					indexBest = i;
					//cout << "Entra: " << all[indexBest]->dist << endl;
				}
			} else {
				if (all[i]->dist >= D){
					if ((betterInFit) || (eqInFit && betterInDist)){
						indexBest = i;
					}
				}
			}
		}
		//Insert best option
		//cout << "Elegido tiene distancia: " << all[indexBest]->dist << endl;
		population.push_back(all[indexBest]);
		acceptedFitness.insert(all[indexBest]->ind.fitness);
		all[indexBest] = all.back();
		all.pop_back();
	}
	//Release memory
	for (int i = 0; i < all.size(); i++){
		delete(all[i]);
	}
}

void MA::initDI(){
	double meanDistance = 0;
	for (int i = 0; i < population.size(); i++){
		for (int j = i + 1; j < population.size(); j++){
			meanDistance += population[i]->ind.getDistance(population[j]->ind);
			//cout << "Distancia: " << population[i]->ind.getDistance(population[j]->ind) << endl;
		}
	}
	meanDistance /= (population.size() * (population.size() - 1)) / 2;
	DI = meanDistance * 0.5;//TODO: Check
}

void MA::run(){
	initPopulation();
	initDI();
	generation = 0;
	while(true){//Infinitas generaciones
		/*int minDistance = INT_MAX;
		int maxDistance = 0;
		cout << "Fitness" << endl;
		for (int i = 0; i < population.size(); i++){
			cout << population[i]->ind.fitness << endl;
			for (int j = i + 1; j < population.size(); j++){
				int d = population[i]->ind.getDistance(population[j]->ind);
				maxDistance = max(maxDistance, d);
				minDistance = min(minDistance, d);
			}
		}
		cout << "Min Distancia: " << minDistance << endl;
		cout << "Max Distancia: " << maxDistance << endl;*/

		//cout << "Generacion " << generation << endl;
		selectParents();
		crossover();
		mutation();
		localSearch();
		replacement();
		generation++;
	}
}

