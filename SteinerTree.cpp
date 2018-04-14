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
	used.insert(u);
	used.insert(v);
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

long long SteinerTree::calculateFitnessComplete(){
	long long F = 0;
	dsu = DSU(SteinerTreeproblem->n);
	FOREACH(e, SteinerTreeproblem->edges) if(I[e->u] && I[e->v]){
		if(!dsu.sameSet(e->u, e->v)){
			dsu.join(e->u, e->v);
			F += e->w;
		}
	}
	dsu.reset();
	return F;
}

bool SteinerTree::calculateFitness(){
	fitness = 0;
	unordered_map<int, int> cnt;
	
	unordered_map<int, vector<pair<int, long long> > > st;
	
	int c = 0;
	FOREACH(e, edges){
		if(!dsu.sameSet(e->u, e->v)){
			dsu.join(e->u, e->v);
			cnt[e->u]++, cnt[e->v]++;
			fitness += e->w;
			st[e->u].push_back(make_pair(e->v, e->w));
			st[e->v].push_back(make_pair(e->u, e->w));
			c++;
		}
	}
	bool connected = true;
	int Su = dsu.used.empty() ? -1 : dsu.getS(*(dsu.used.begin()));
	FOREACH(v, dsu.used) if(dsu.getS(*v) != Su) connected = false;
	int t = 0;
	for(int u = 0; u < (SteinerTreeproblem->n); u++) if(I[u]) t++;
	connected = (c == t - 1);
	/*if(connected && c != t - 1){
		printf("NO = %d, %d\n", c, t);
		getchar();
	}*/
	if(connected){
		queue<int> q;
		FOREACH(v, dsu.used) if(cnt[*v] == 1 && !(SteinerTreeproblem->fixed)[*v]) q.push(*v);
		while(!q.empty()){
			int u = q.front(); q.pop();
			I[u] = false;
			FOREACH(v, (SteinerTreeproblem->adj)[u]) if(I[v->first]) edges.erase(edge(min(u, v->first), max(u, v->first), v->second));
			FOREACH(v, st[u]) if(I[v->first]){
				cnt[v->first]--;
				fitness -= v->second;
				if(cnt[v->first] == 1 && !(SteinerTreeproblem->fixed)[v->first]) q.push(v->first);
			}
		}
	}
	else fitness = (long long)1e15;
	/*if(fitness == 0){
		for(int u = 0; u < (SteinerTreeproblem->n); u++) if((SteinerTreeproblem->fixed)[u]) printf("%d ", I[u] ? 1 : 0);
		printf("\n");
		printf("connected = %d\n", connected ? 1 : 0);
		FOREACH(e, edges) printf("%d %d\n", e->u, e->v);
		getchar();
	}*/
	dsu.reset();
	return connected;
}

void SteinerTree::insert(int u){
	if(!I[u]){
		I[u] = true;
		FOREACH(vw, (SteinerTreeproblem->adj)[u]) if(I[vw->first]){
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
		FOREACH(vw, (SteinerTreeproblem->adj)[u]) if(I[vw->first]){
			int v = vw->first;
			long long w = vw->second;
			edges.erase(edge(min(u, v), max(u, v), w));
		}
		calculateFitness();
	}
}

long long best = 1e16;
void SteinerTree::hillClimbing(){
	vector<int> p(SteinerTreeproblem->n);
	for(int i = 0; i < SteinerTreeproblem->n; i++) p[i] = i;
	calculateFitness();
	int noImprove = 0;
	while(true){
		random_shuffle(p.begin(), p.end());
		bool done = false;
		for(int i = 0; i < SteinerTreeproblem->n; i++){
			if((SteinerTreeproblem->fixed)[p[i]]) continue;
			best = min(best, fitness);
			printf("Fitness = %lld %lld\n", fitness, best);
			long long F = fitness;
			vector<bool> nI = I;
			if(I[p[i]]){
				erase(p[i]);
				if(fitness < F){
					done = false;
					continue;
				}
			}
			else{
				insert(p[i]);
				if(fitness < F){
					done = false;
					continue;
				}
			}
			noImprove++;
			reset(nI);
			printf("noImprove = %d\n", noImprove);
			if(noImprove > 100) return;
		}
		if(done) break;
	}
}

void SteinerTree::localSearch(){
	while(true){
		hillClimbing();
		for(int i = 0; i < 20; i++){
			int u = rand()%(SteinerTreeproblem->n);
			FOREACH(v, (SteinerTreeproblem->adj)[u]) if(I[v->first]){
				insert(u);
				break;
			}
		}
	}
}

void SteinerTree::reset(vector<bool> &nI){
	I = nI;
	edges.clear();
	for(int u = 0; u < (SteinerTreeproblem->n); u++)
		if(nI[u]) FOREACH(v, (SteinerTreeproblem->adj)[u]) if(v->first < u && I[v->first]) edges.insert(edge(v->first, u, v->second));
	calculateFitness();
}

SteinerTreeProblem::SteinerTreeProblem(const string &fileName){
	ifstream ifs;
	ifs.open(fileName.c_str());
	string s;
	ifs >> s;
	ifs >> s;
	ifs >> s;
	ifs >> n;
	ifs >> s;
	ifs >> m;
	int u, v;
	long long w;
	adj.resize(n);
	fixed.resize(n);
	for(int i = 0; i < m; i++){
		ifs >> s >> u >> v >> w; u--, v--;
		edges.push_back(edge(u, v, w));
		adj[u].push_back(make_pair(v, w));
		adj[v].push_back(make_pair(u, w));
	}
	sort(edges.begin(), edges.end());
	ifs >> s;
	ifs >> s;
	ifs >> s;
	ifs >> s;
	int T; ifs >> T;
	for(int i = 0; i < T; i++){
		ifs >> s >> u; u--;
		fixed[u] = true;
	}
}

void SteinerTree::restart(){
	I.resize(SteinerTreeproblem->n);
	for(int i = 0; i < (int)I.size(); i++) I[i] = true;
	for(int u = 0; u < SteinerTreeproblem->n; u++)
		FOREACH(v, (SteinerTreeproblem->adj)[u]) if(v->first < u){
			edges.insert(edge(v->first, u, v->second));
		}
	dsu = DSU((int)I.size());
	localSearch();
}

