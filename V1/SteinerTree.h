#ifndef __STEINER_TREE_H__
#define __STEINER_TREE_H__

#include <bits/stdc++.h>
using namespace std;
#define FOREACH(i, v) for (__typeof((v).begin()) i = (v).begin(); i != (v).end(); i++)

extern volatile bool finished;

struct edge{
	long long w;
	int u, v;
	edge(int u, int v, long long w){
		this->u = u;
		this->v = v;
		this->w = w;
	}
	edge(){
	
	}
};

void printBest();

class SteinerTreeProblem{
	public:
		SteinerTreeProblem();
		~SteinerTreeProblem(){
		
		}
		vector<vector<pair<int, long long> > > adj;
		vector<edge> edges;
		vector<bool> fixed;
		vector<int> fs;
		int n, m;

};


class DSU{
	public:
		DSU(int n);
		DSU(){
		
		}
		~DSU(){
		
		}
		int getS(int u);
		bool sameSet(int u, int v);
		void join(int u, int v);
		void reset();
		vector<int> p;
		unordered_set<int> used;
};


class SteinerTree{
	public:
		SteinerTree(){
		
		}
		~SteinerTree(){
		
		}
		void restart();
		void reset(vector<bool> &nI);
		long long calculateFitnessComplete();
		int getDistance(SteinerTree &ind);
		//ostream& operator<< (ostream &os, const SteinerTree &ST);
		//void print(ostream &os) const;
		void getMST(unordered_map<int, vector<pair<int, long long> > > &mst);
		void dependentMutation(double pm);
		void uniformMutation(double pm);
		void pathMutation(int k);
		void localSearch();
		void hillClimbing();
		void dependentCrossover(SteinerTree &ind);
		void uniformCrossover(SteinerTree &ind);
		void addCrossover(SteinerTree &ind2);
		bool calculateFitness();
		void insert(int u);
		void erase(int u);
		void print(unordered_map<int, vector<pair<int, long long> > > &mst);
		bool isCorrect(unordered_map<int, vector<pair<int, long long> > > &mst);

		vector<bool> I;
		long long fitness;
		set<edge> edges;
		static SteinerTreeProblem *SteinerTreeproblem;


};

#endif
