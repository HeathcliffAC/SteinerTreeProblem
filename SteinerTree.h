#include <bits/stdc++.h>
using namespace std;
#define FOREACH(i, v) for (__typeof((v).begin()) i = (v).begin(); i != (v).end(); i++)

struct edge{
	long long w;
	int u, int v;
	edge(int u, int v, long long w){
		this->u = u;
		this->v = v;
		this->w = w;
	}
	edge(){
	
	}
};

friend bool operator< (const edge& a, const edge& b){
	return a.w < b.w;
}

class SteinerTreeProblem{
	public:
		SteinerTreeProblem(const string &fileName);
		SteinerTreeProblem(){
		
		}
		~SteinerTreeProblem(){
		
		}
		vector<vector<pair<int, long long> > > adj;
		vector<edge> edges;
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
		unordered_set used;
};

DSU dsu; 

class SteinerTree{
	public:
		SteinerTree(){
		
		{
		~SteinerTree(){
		
		}
		void restart();
		double getDistance(SteinerTree &ind);
		friend ostream& operator<< (ostream &os, const SteinerTree &ST);
		void print(ostream &os) const;
		void dependentMutation(double pm);
		void localSearch();
		void dependentCrossover(SteinerTree *ind);
		void calculateFitness();
		void insert(int u);
		void erase(int u);
		vector<bool> I;
		long long fitness;
		set<edge> edges;
		SteinerTreeProblem *SteinerTreeproblem;
};
