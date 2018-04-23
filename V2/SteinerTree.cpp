#include <signal.h>
#include "SteinerTree.h"
#include "utils.h"
using namespace std;

bool volatile finished = false;
SteinerTreeProblem* SteinerTree::SteinerTreeproblem;
DSU dsu; 
bool operator< (const edge& a, const edge& b){
	if(a.w < b.w) return true;
	else if(a.w > b.w) return false;
	else return make_pair(a.u, a.v) < make_pair(b.u, b.v);
}

///////////////////// Funciones de un union find disjoints sets ///////////////////////////////
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


//////////////////////////// Funciones del individuo ////////////////////////////////

// Calcula el fitness "por fuerza bruta"
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

bool continuar = false;

// Calcula fitness, devuelve si esta conectado, poda las hojas no necesarias
bool SteinerTree::calculateFitness(){
	evaluateMinDistances();
	if(!continuar) return true;
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
	int Su = -1;
	Su = dsu.getS((SteinerTreeproblem->fs)[0]);
	FOREACH(u, SteinerTreeproblem->fs) if(dsu.getS(*u) != Su) connected = false;
	//FOREACH(v, dsu.used) if(dsu.getS(*v) != Su) connected = false;
	//int t = 0;
	//for(int u = 0; u < (SteinerTreeproblem->n); u++) if((SteinerTreeproblem->fixed)[u] && dsu.getS(u) != Su) connected = false;
	//connected = (c == t - 1);
	if(connected){
		queue<int> q;
		//for(int u = 0; u < (SteinerTreeproblem->n); u++) if(cnt.count(u) == 0 && I[u]) I[u] = false;
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
	dsu.reset();
	return connected;
}

// Inserta el nodo u al individuo
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

void SteinerTree::insert(vector<int> us){
	FOREACH(u, us){
		if(!I[*u]){
			I[*u] = true;
			FOREACH(vw, (SteinerTreeproblem->adj)[*u]) if(I[vw->first]){
				int v = vw->first;
				long long w = vw->second;
				edges.insert(edge(min(*u, v), max(*u, v), w));
			}
		}
		calculateFitness();
	}
}


// Elimina el nodo u del individuo
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

long long Globalbest = 1e16;
SteinerTree bestI;

void printBest(){
	unordered_map<int, vector<pair<int, long long> > > mst;
	bestI.getMST(mst);
	//bestI.print(mst);
	bool correct = bestI.isCorrect(mst);
	//printf("%s\n", correct ? "YES" : "NO");
	if(!correct){
		printf("NO\n");
		//getchar();
	}
	else bestI.print(mst);
}



// Primera implementacion de una busqueda por escalada
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
			if (fitness < Globalbest){
				sigset_t sign;
				sigemptyset(&sign);
				sigaddset(&sign, SIGTERM);
				sigprocmask(SIG_BLOCK, &sign, NULL);
				//printf("Fitness = %lld\n", fitness);
				Globalbest = fitness;
				bestI = *this;
				sigprocmask(SIG_UNBLOCK, &sign, NULL);
				//printBest();
			}
			//best = min(best, fitness);
			//printf("Fitness = %lld %lld\n", fitness, best);
			long long F = fitness;
			vector<bool> nI = I;
			if(I[p[i]]){
				erase(p[i]);
				if(fitness < F){
					noImprove = 0;
					done = false;
					continue;
				}
			}
			else{
				/*vector<int> vc; vc.push_back(p[i]);
				for(int k = 0; k < 2; k++) vc.push_back(rand()%((int)p.size()));
				insert(vc);
				*/
				
				insert(p[i]);
				if(fitness < F){
					noImprove = 0;
					done = false;
					continue;
				}
			}
			//if(rand()%2 == 0 && I[p[i]]) continue;
			noImprove++;
			reset(nI);
			if(noImprove > 5) return;
		}
	}
}

bool SteinerTree::isCorrect(unordered_map<int, vector<pair<int, long long> > > &mst){
	int cnt = 0;
	vector<bool> visited(I.size());
	bool ans = true;
	long long fitnessAux = 0;
	FOREACH(u, SteinerTreeproblem->fs){
		if(!I[*u]){
			printf("u = %d\n", *u);
			ans = false;
		}
		if(!visited[*u]){
			int totalNodes = 0, totalEdges = 0;
			queue<int> q; q.push(*u);
			visited[*u] = true;
			while(!q.empty()){
				int v = q.front(); q.pop();
				totalNodes++;
				FOREACH(w, mst[v]){
					if(!visited[w->first]){
						q.push(w->first);
						visited[w->first] = true;
						fitnessAux += w->second;
					}
					if(w->first < v) totalEdges++;
				}
			}
			if(totalEdges != totalNodes - 1){
				printf("t = %d, t = %d\n", totalEdges, totalNodes);
				ans = false;
			}
			cnt++;
		}
	}
	if(fitnessAux > Globalbest){
		printf("fitnessAux = %lld, fitness = %lld\n", fitnessAux, Globalbest);
		ans = false;
	}
	else Globalbest = fitnessAux;
	if(cnt != 1){
		printf("cnt = %d\n", cnt);
		ans = false;
	}
	return ans;
}

void SteinerTree::localSearch(){
	return;
	hillClimbing();
	printf("fitness = %lld %lld\n", fitness, Globalbest);
}

void SteinerTree::reset(vector<bool> &nI){
	I = nI;
	edges.clear();
	for(int u = 0; u < (SteinerTreeproblem->n); u++)
		if(nI[u]) FOREACH(v, (SteinerTreeproblem->adj)[u]) if(v->first < u && I[v->first]) edges.insert(edge(v->first, u, v->second));
	calculateFitness();
}

SteinerTreeProblem::SteinerTreeProblem(){
	//ifstream ifs;
	//ifs.open(fileName.c_str());
	string s;
	cin >> s;
	cin >> s;
	cin >> s;
	cin >> n;
	cin >> s;
	cin >> m;
	int u, v;
	long long w;
	adj.resize(n);
	fixed.resize(n);
	for(int i = 0; i < m; i++){
		cin >> s >> u >> v >> w; u--, v--;
		edges.push_back(edge(u, v, w));
		adj[u].push_back(make_pair(v, w));
		adj[v].push_back(make_pair(u, w));
	}
	sort(edges.begin(), edges.end());
	cin >> s;
	cin >> s;
	cin >> s;
	cin >> s;
	int T; cin >> T;
	for(int i = 0; i < T; i++){
		cin >> s >> u; u--;
		fixed[u] = true;
		fs.push_back(u);
	}
}

vector<int> p, sz;
int mx;

void init(vector<bool> &I){
	p.resize(I.size());
	sz.resize(I.size());
	int len = p.size();
	for(int i = 0; i < len; i++) p[i] = i;
	for(int i = 0; i < len; i++) sz[i] = (I[i] ? 1 : 0);
}

int getS(int u){
	if(p[u] == u) return u;
	else return p[u] = getS(p[u]);
}

bool sameS(int u, int v){
	return getS(u) == getS(v);
}

void join(int u, int v){
	u = getS(u), v = getS(v);
	if(u == v) return;
	p[u] = v;
	sz[v] += sz[u];
	mx = max(mx, sz[v]);
}



bool first = true;
void SteinerTree::restart(){
	I.resize(SteinerTreeproblem->n, false);
	for(int i = 0; i < (int)(SteinerTreeproblem->fs).size(); i++)
		I[(SteinerTreeproblem->fs)[i]] = true;
	for(int i = 0; i < (SteinerTreeproblem->n); i++) if(rand()/(RAND_MAX + 1.0) < 0.15) I[i] = true;
	dsu = DSU((int)I.size());
	calculateFitness();
	return;
	
	init(I);
	int ttl = (int)(SteinerTreeproblem->fs).size();
	mx = 1;
	
	for(int i = 0; i < (int)(SteinerTreeproblem->fs).size(); i++){
		int u = (SteinerTreeproblem->fs)[i];
		FOREACH(v, (SteinerTreeproblem->adj)[u]) if(I[v->first]){
			edges.insert(edge(min(u, v->first), max(u, v->first), v->second));
			join(u, v->first);
		}
	}
	
	vector<int> vc;
	for(int i = 0; i < (SteinerTreeproblem->n); i++) if(!I[i]) vc.push_back(i);
	
	while(mx != ttl){
		int id = rand()%((int)vc.size());
		int u = vc[id];
		I[u] = true;
		FOREACH(v, (SteinerTreeproblem->adj)[u]) if(I[v->first]){
			edges.insert(edge(min(u, v->first), max(u, v->first), v->second));
			join(u, v->first);
		}
		swap(vc[id], vc[(int)vc.size() - 1]);
		vc.pop_back();
	}
	dsu = DSU((int)I.size());
	
	
	
	localSearch();
}

//Se incluyen algunos del otro y se dejan todos los que estan
/*void SteinerTree::crossover(SteinerTree &ind2){
	for (int i = 0; i < I.size(); i++){
		if (I[i]){
			if (generateRandomDouble0_Max(1) < 0.5){
				ind2.I[i] = true;
			}
		}
	}
	for (int i = 0; i < I.size(); i++){
		if (ind2.I[i]){
			if (generateRandomDouble0_Max(1) < 0.5){
				I[i] = true;
			}
		}
	}
}*/

/*
void SteinerTree::mutate(double pm){
	for (int i = 0; i < I.size(); i++){
		if (generateRandomDouble0_Max(1) < pm){
			I[i] = true;
		}
	}
}*/

void SteinerTree::addCrossover(SteinerTree &ind2){
	for (int i = 0; i < I.size(); i++){
		if (I[i]){
			if (generateRandomDouble0_Max(1) < 0.5){
				ind2.I[i] = true;
			}
		}
	}
	for (int i = 0; i < I.size(); i++){
		if (ind2.I[i]){
			if (generateRandomDouble0_Max(1) < 0.5){
				I[i] = true;
			}
		}
	}
	reset(I);
	ind2.reset(ind2.I);
}

void SteinerTree::uniformCrossover(SteinerTree &ind){
	for(int i = 0; i < (SteinerTreeproblem->n); i++)
		if(rand()%2 == 0) swap(I[i], (ind.I)[i]);
	
	//reset(I);
	//ind.reset(ind.I);
	calculateFitness();
	ind.calculateFitness();
}

int SteinerTree::getDistance(SteinerTree &ind){
	int ans = 0;
	for(int i = 0; i < (SteinerTreeproblem->n); i++)
		if(I[i] != (ind.I)[i]) ans++;
	return ans;
}

void SteinerTree::pathMutation(int k){
	int u = 1;
	while(!I[u]) u = rand()%(SteinerTreeproblem->n);
	for(int i = 0; i < k; i++){
		I[u] = true;
		int id = rand()%((int)(SteinerTreeproblem->adj)[u].size());
		u = (SteinerTreeproblem->adj)[u][id].first;
	}
	reset(I);
}

void SteinerTree::uniformMutation(double pm){
	long long prevFitness = fitness;
	bool changed = false;
	vector<bool> prevI = I;
	for(int i = 0; i < (SteinerTreeproblem->n); i++) if(!(SteinerTreeproblem->fixed)[i])
		if(rand()/(RAND_MAX + 1.0) < pm){
			changed = true;
			I[i] = !I[i];
		}
	//reset(I);
	if (changed){
		calculateFitness();
		//cout << "Comp: " << prevFitness << " " << fitness << endl;
		if (prevFitness < fitness){
			I = prevI;
			fitness = prevFitness;
		}
	}
}

void SteinerTree::dependentMutation(double pm){
	//cout << "empieza mut" << endl;
	for (int i = 0; i < 10; i++){
		uniformMutation(1.0 / (SteinerTreeproblem->n));
	}
	//cout << "termina" << endl;
}

void SteinerTree::dependentCrossover(SteinerTree &ind){
	uniformCrossover(ind);
}

// Guarda en mst el minimum spanning tree del individuo
void SteinerTree::getMST(unordered_map<int, vector<pair<int, long long> > > &mst){
	edges.clear();
	for(int u = 0; u < (SteinerTreeproblem->n); u++) if(I[u]){
		FOREACH(v, (SteinerTreeproblem->adj)[u]) if(u < v->first && I[v->first]) edges.insert(edge(u, v->first, v->second));
	}
	mst.clear();
	FOREACH(e, edges){
		if(!dsu.sameSet(e->u, e->v)){
			dsu.join(e->u, e->v);
			mst[e->u].push_back(make_pair(e->v, e->w));
			mst[e->v].push_back(make_pair(e->u, e->w));
		}
	}
	dsu.reset();
}

void SteinerTree::print(unordered_map<int, vector<pair<int, long long> > > &mst){
    unordered_map<int, vector<pair<int, long long> > >::iterator it;

    long long cost = 0;
    for(it = mst.begin(); it != mst.end(); it ++)
    	FOREACH(v, it->second) cost += (v->second);
    
    cost = cost/2;

    cout << "VALUE " << cost << endl;

    for(it = mst.begin(); it != mst.end(); it ++)
    	FOREACH(v, it->second)
        if(it->first < v->first)
            cout << it->first + 1 << " " << v->first + 1 << endl;
}

void SteinerTreeProblem::dijkstra(vector<int> &u, vector<long long> &dist, vector<int> &p, vector<bool> &tree, priority_queue<pair<long long, int>, vector<pair<long long, int> >, greater<pair<long long, int> > > &pq, vector<bool> &I){
	if(p.empty()){
		long long INF = 1e15;
		p = vector<int>(n);
		dist = vector<long long>(n, INF);
		tree = vector<bool>(n, false);
	}
	
	//priority_queue<pair<long long, int>, vector<pair<long long, int> >, greater<pair<long long, int> > > pq;
	for (int i = 0; i < u.size(); i++){
		tree[u[i]] = true;
		dist[u[i]] = 0; pq.push(make_pair(dist[u[i]], u[i]));
	}
	while(!pq.empty()){
		pair<long long, int> tmp = pq.top(); pq.pop();
		long long d = tmp.first;
		int n = tmp.second;
		if(dist[n] != d) continue;
		FOREACH(vw, adj[n]){
			int v = vw->first;
			long long w = vw->second;
			if(!tree[v] && dist[v] > dist[n] + w){
				dist[v] = dist[n] + w;
				p[v] = n;
				pq.push(make_pair(dist[v], v));
			}
		}
		if (I[n] && (!tree[n])){
			return;
		}
	}
}

long long fitnessBest = 1e15;

void SteinerTree::evaluateMinDistances(){
	long long best = 1e15;
	vector<bool> Ibest;
	int cnt_break = 0;
	vector<bool> bestTree;
	while(true){
		int u = (SteinerTreeproblem->fs)[rand()%((int)(SteinerTreeproblem->fs).size())];
		fitness = 0;
		vector<long long> dist;
		vector<int> p;
		vector<bool> tree;
		vector<int> nodes;
		nodes.push_back(u);
		priority_queue<pair<long long, int>, vector<pair<long long, int> >, greater<pair<long long, int> > > pq;
		SteinerTreeproblem->dijkstra(nodes, dist, p, tree, pq, I);
		vector<int> cnt(SteinerTreeproblem->n, 0);
		while(true){
			long long mn = 1e15;
			u = -1;
			for(int i = 0; i < (SteinerTreeproblem->n); i++){
				if(I[i] && !tree[i] && dist[i] < mn) mn = dist[i], u = i;
			}
			if(u == -1) break;
			fitness += mn;
			vector<int> path;
			while(!tree[u]){
				path.push_back(u);
				cnt[u]++, cnt[p[u]]++;
				u = p[u];
			}
			reverse(path.begin(), path.end());
			//for(int i = 0; i < (int)path.size(); i++) SteinerTreeproblem->dijkstra(path[i], dist, p, tree);
			SteinerTreeproblem->dijkstra(path, dist, p, tree, pq, I);
			if (finished){
				printBest();
				exit(0);
			}
		}
		if(best > fitness){
			best = fitness;
			Ibest = I;
			bestTree = tree;
			cnt_break = 0;
			if(fitness < Globalbest){
				//cout << "Va por " << fitness << endl;
				sigset_t sign;
				sigemptyset(&sign);
				sigaddset(&sign, SIGTERM);
				//sigaddset(&sign, SIGALRM);
				sigprocmask(SIG_BLOCK, &sign, NULL);
				Globalbest = fitness;
				bestI.I = bestTree;
				sigprocmask(SIG_UNBLOCK, &sign, NULL);
			}
		}
		else cnt_break++;
		if(cnt_break > 3) break;
		for(int u = 0; u < (SteinerTreeproblem->n); u++){
			if((SteinerTreeproblem->fixed)[u] || cnt[u] >= 3) I[u] = true; 
			else I[u] = false;
		}
	}
	I = Ibest;
	fitness = best;
	/*if(fitness < Globalbest){
		Globalbest = fitness;
		bestI.I = bestTree;
	}*/
	//printf("best = %lld, fitness = %lld\n", Globalbest, fitness);
}


