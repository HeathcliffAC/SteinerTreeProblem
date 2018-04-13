#include "SteinerTree.h"
using namespace std;
DSU::DSU(int n){
	p.resize(n + 1);
	for(int i = 0; i <= n; i++) p[i] = i;
}

int DSU::getS(int u){
	return (p[u] == u ? u : (p[u] = getS(p[u])));
}

bool DSU::sameSet(int u, int v){
	return getS(u) == getS(v);
}

void DSU::join(int u, int v){
	used.insert(u);
	used.insert(v);
	u = getS(u);
	v = getS(v);
	if(u != v) p[u] = v;
}

void DSU::reset(){
	FOREACH(u, used) p[*u] = *u;
	used.clear();
}

bool SteinerTree::calculateFitness(){
	fitness = 0;
	unordered_map<int, int> cnt;
	FOREACH(e, edges){
		if(!dsu.sameSet(e->u, e->v)){
			dsu.join(e->u, e->v);
			cnt[e->u]++, cnt[e->v]++;
			fitness += e->w;
		}
	}
	bool connected = true;
	int Su = dsu.getS(*dsu.used.begin());
	FOREACH(v, dsu.used) if(dsu.getS(*v) != Su) connected = false;
	dsu.reset();
	if(connected){
		queue<int> q;
		FOREACH(v, dsu.used) if(cnt[*v] == 1 && !(SteinerTreeproblem->fixed)[*v]) q.push(*v);
		while(!q.empty()){
			int u = q.front();
			FOREACH(v, (SteinerTreeproblem->adj)[u]) if(I[v->first]){
				cnt[v->first]--;
				fitness -= v->second;
				if(cnt[v->first] == 1 && !(SteinerTreeproblem->fixed)[v->first]) q.push(v->first);
			}
		}
	}
	return connected;
}

void SteinerTree::insert(int u){
	if(!I[u]){
		I[u] = true;
		FOREACH(vw, SteinerTreeproblem->adj[u]) if(I[vw->first]){
			int v = vw->first;
			long long w = vw->second;
			edges.insert(edge(min(u, v), max(u, v), w));
		}
		calculateFitness();
	}
}

void SteinerTree::erase(int u){
	if(I[u]){
		I[u] = false;
		FOREACH(vw, SteinerTreeproblem->adj[u]) if(I[vw->first]){
			int v = vw->first;
			long long w = vw->second;
			edges.erase(edge(min(u, v), max(u, v), w));
		}
		calculateFitness();
	}
}

void SteinerTree::hillClimbing(){
	vector<int> p(SteinerTreeproblem->n);
	for(int i = 0; i < SteinerTreeproblem->n; i++) p[i] = i;
	calculateFitness();
	while(true){
		random_shuffle(p.begin(), p.end());
		for(int i = 0; i < SteinerTreeproblem->n; i++){
			printf("Fitness = %lld\n", fitness);
			long long F = fitness;
			if(I[p[i]]){
				erase(p[i]);
				if(fitness < F) continue;
				insert(p[i]);
			}
			else{
				insert(p[i]);
				if(fitness < F) continue;
				erase(p[i]);
			}
		}
	}
}

void SteinerTree::localSearch(){
	hillClimbing();
}
