#include "SteinerTree.h"

DSU::DSU(int n){
	p.resize(n + 1);
	for(int i = 0; i <= n; i++) p[i] = i;
}

int DSU::getS(int u){
	return p[u] == u ? u : (p[u] = getS(p[u]));
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

void SteinerTree::calculateFitness(){
	fitness = 0;
	FOREACH(e, edges){
		if(!dsu.sameSet(e->u, e->v)){
			dsu.join(e->u, e->v);
			fitness += e->w;
		}
	}
	dsu.reset();
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


