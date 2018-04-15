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

bool operator< (const edge& a, const edge& b){
	if(a.w < b.w) return true;
	else if(a.w > b.w) return false;
	else return make_pair(a.u, a.v) < make_pair(b.u, b.v);
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

DSU dsu; 

class SteinerTree{
	public:
		SteinerTree(){
		
		}
		~SteinerTree(){
		
		}
		void restart();
		void reset(vector<bool> &nI);
		long long calculateFitnessComplete();
		double getDistance(SteinerTree &ind);
		//ostream& operator<< (ostream &os, const SteinerTree &ST);
		//void print(ostream &os) const;
		void getMST(unordered_map<int, vector<pair<int, long long> > > &mst);
		void dependentMutation(double pm);
		void uniformMutation(double pm);
		void pathMutation(int k);
		void localSearch();
		void hillClimbing();
		void dependentCrossover(SteinerTree* ind);
		void uniformCrossover(SteinerTree* ind);
		bool calculateFitness();
		void insert(int u);
		void erase(int u);
		vector<bool> I;
		long long fitness;
		set<edge> edges;
		SteinerTreeProblem *SteinerTreeproblem;
};
