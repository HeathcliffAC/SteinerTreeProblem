#ifndef __STEINER_TREE_H__
#define __STEINER_TREE_H__

#include <bits/stdc++.h>
using namespace std;
#define FOREACH(i, v) for (__typeof((v).begin()) i = (v).begin(); i != (v).end(); i++)

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


class SteinerTreeProblem{
	public:
		SteinerTreeProblem(const string &fileName);
		SteinerTreeProblem(){
		
		}
		~SteinerTreeProblem(){
		
		}
		vector<vector<pair<int, long long> > > adj;
		vector<edge> edges;
		vector<bool> fixed;
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
		//double getDistance(SteinerTree &ind);
		//ostream& operator<< (ostream &os, const SteinerTree &ST);
		//void print(ostream &os) const;
		//void dependentMutation(double pm);
		void localSearch();
		void crossover(SteinerTree &ind2);//Modifies both "this" and "ind2"
		void mutate(double pm);
		void hillClimbing();
		int getDistance(SteinerTree &ind2);
		//void dependentCrossover(SteinerTree *ind);
		bool calculateFitness();
		void insert(int u);
		void erase(int u);
		vector<bool> I;
		long long fitness;
		set<edge> edges;
		static SteinerTreeProblem *SteinerTreeproblem;
};

#endif
