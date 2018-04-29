int getRandomInteger0_N(int n);//get random integer [0, n]
double generateRandomDouble0_Max(double maxValue);
#include <stdlib.h>

int getRandomInteger0_N(int n){ 
  return (int) ((n + 1.0)*rand()/(RAND_MAX+1.0));
}

double generateRandomDouble0_Max(double maxValue){
	return (double)(rand()) / RAND_MAX * maxValue;
}
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
		void dijkstra(vector<int> &u, vector<long long> &dist, vector<int> &p, vector<bool> &done, priority_queue<pair<long long, int>, vector<pair<long long, int> >, greater<pair<long long, int> > > &pq, vector<bool> &I);
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
		void reset();
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
		bool calculateFitness2();
		void evaluateMinDistances();
		void insert(int u);
		void insert(vector<int> us);
		void erase(int u);
		void print(unordered_map<int, vector<pair<int, long long> > > &mst);
		bool isCorrect(unordered_map<int, vector<pair<int, long long> > > &mst);
		vector<bool> I, IComplete;
		long long fitness;
		vector<edge> edges;
		static SteinerTreeProblem *SteinerTreeproblem;
};

#endif
#include <signal.h>
using namespace std;

int currentTol;
long long Globalbest = 1e16;
SteinerTree bestI;
int generation = 0;
bool volatile finished = false;
SteinerTreeProblem* SteinerTree::SteinerTreeproblem;
DSU dsu; 
long long maxThreshold;
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

bool SteinerTree::calculateFitness2(){
	fitness = 0;
	int cnt[SteinerTreeproblem->n];
	bzero(cnt, sizeof(cnt));
	vector<pair<int, long long> > st[SteinerTreeproblem->n];
	
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
	queue<int> q;
	FOREACH(v, dsu.used) if(cnt[*v] == 1 && !(SteinerTreeproblem->fixed)[*v]) q.push(*v);
	while(!q.empty()){
		int u = q.front(); q.pop();
		I[u] = false;
		cnt[u] = 0;
		FOREACH(v, st[u]) if(I[v->first]){
			cnt[v->first]--;
			fitness -= v->second;
			if(cnt[v->first] == 1 && !(SteinerTreeproblem->fixed)[v->first]) q.push(v->first);
		}
	}
	dsu.reset();

	if(fitness < Globalbest){
		sigset_t sign;
		sigemptyset(&sign);
		sigaddset(&sign, SIGTERM);
		//sigaddset(&sign, SIGALRM);
		sigprocmask(SIG_BLOCK, &sign, NULL);
		bestI.I.clear();
		for(int u = 0; u < (SteinerTreeproblem->n); u++){
			if((SteinerTreeproblem->fixed)[u] || cnt[u] >= 1){
				bestI.I.push_back(true);
			} else {
				bestI.I.push_back(false);
			}
		}
		Globalbest = fitness;
		sigprocmask(SIG_UNBLOCK, &sign, NULL);
	}
	return true;
}

bool SteinerTree::calculateFitness(){
	evaluateMinDistances();
	return true;
}

// Inserta el nodo u al individuo
void SteinerTree::insert(int u){
}

void SteinerTree::insert(vector<int> us){
}


// Elimina el nodo u del individuo
void SteinerTree::erase(int u){
}

void printBest(){
	unordered_map<int, vector<pair<int, long long> > > mst;
	bestI.getMST(mst);
	//bestI.print(mst);
	bool correct = bestI.isCorrect(mst);
	//printf("%s\n", correct ? "YES" : "NO");
	if(!correct){
		printf("NO\n");
	}
	else bestI.print(mst);
}



// Primera implementacion de una busqueda por escalada
void SteinerTree::hillClimbing(){
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

/*
void SteinerTree::reset(vector<bool> &nI){
	I = nI;
	edges.clear();
	for(int u = 0; u < (SteinerTreeproblem->n); u++)
		if(nI[u]) FOREACH(v, (SteinerTreeproblem->adj)[u]) if(v->first < u && I[v->first]) edges.insert(edge(v->first, u, v->second));
	calculateFitness2();
}*/

void SteinerTree::reset(){
	edges.clear();
	for(int u = 0; u < (SteinerTreeproblem->n); u++)
		if(I[u]) FOREACH(v, (SteinerTreeproblem->adj)[u]) if(v->first < u && I[v->first]) edges.push_back(edge(v->first, u, v->second));
	sort(edges.begin(), edges.end());
	calculateFitness2();
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
		//if(w == 0) while(true) printf("u = %d, v = %d\n", u, v);
	}
	sort(edges.begin(), edges.end());
	maxThreshold = edges[edges.size() / 2].w * 2;
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
	IComplete.resize(SteinerTreeproblem->n, false);
	for(int i = 0; i < (int)(SteinerTreeproblem->fs).size(); i++)
		I[(SteinerTreeproblem->fs)[i]] = true;
	for(int i = 0; i < (SteinerTreeproblem->n); i++) 
		if(rand()/(RAND_MAX + 1.0) < 0.15) I[i] = true;
	dsu = DSU((int)I.size());

	//I[45-1] = I[170-1] = I[10-1] = I[5-1] = true; //I[2-1] = true;	

	calculateFitness();
	return;
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
	//cout << "Dist1: " << ans << endl;
	int ans2 = 0;
	for(int i = 0; i < (SteinerTreeproblem->n); i++)
		if(IComplete[i] != (ind.IComplete)[i]) ans2++;
	//cout << "Dist2: " << ans2 << endl;
	return ans;
}

void SteinerTree::pathMutation(int k){
}

void SteinerTree::uniformMutation(double pm){
	long long prevFitness = fitness;
	bool changed = false;
	vector<bool> prevI = I;
	vector<bool> prevIComplete = IComplete;
	for(int i = 0; i < (SteinerTreeproblem->n); i++) if(!(SteinerTreeproblem->fixed)[i])
		if(rand()/(RAND_MAX + 1.0) < pm){
			changed = true;
			I[i] = !I[i];
		}
	/*vector<int> ones, zeros;
	for(int i = 0; i < (SteinerTreeproblem->n); i++) 
		if(!(SteinerTreeproblem->fixed)[i])
			if (I[i])
				ones.push_back(i);
			else
				zeros.push_back(i);
	if (ones.size() == 0) return;
	if (zeros.size() == 0) return;*/
	/*cout << "Solucion" << endl;
	for (int i = 0; i < ones.size(); i++){
		cout << ones[i] << " ";
	}
	cout << endl;*/
/*	int r1 = zeros[rand() % zeros.size()];
	I[zeros[r1]] = true;
	while(ones.size() > 4){
		int index = rand() % ones.size();
		int r2 = ones[index];
		I[r2] = false;
		ones[index] = ones.back();
		ones.pop_back();
	}
	ones.clear();
	for(int i = 0; i < (SteinerTreeproblem->n); i++) 
		if(!(SteinerTreeproblem->fixed)[i])
			if (I[i])
				ones.push_back(i);
			else
				zeros.push_back(i);*/
	//cout << "Cantidad: " << ones.size() << endl;
	/*int r1 = zeros[rand() % zeros.size()], r2 = ones[rand() % ones.size()];//Warning
	I[r1] = !I[r1];
	I[r2] = !I[r2];*/
	//changed = true;
	//reset(I);
	if (changed){
		calculateFitness();
		//cout << "Comp: " << prevFitness << " " << fitness << endl;
		if (prevFitness < fitness){
			I = prevI;
			IComplete = prevIComplete;
			fitness = prevFitness;
		}
	}
}

void SteinerTree::dependentMutation(double pm){
	//cout << "empieza mut" << endl;
	int limit = (generation < 100)?0:5;
	//limit = 1;
	for (int i = 0; i < limit; i++){
		if (SteinerTreeproblem->n - SteinerTreeproblem->fs.size() != 0){
			uniformMutation(2.0 / (SteinerTreeproblem->n - SteinerTreeproblem->fs.size()));
		}
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
		FOREACH(v, (SteinerTreeproblem->adj)[u]) if(u < v->first && I[v->first]) edges.push_back(edge(u, v->first, v->second));
	}
	sort(edges.begin(), edges.end());
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

vector<int> hsh_r;

long long counter = 0;

set<pair<long long, int> > order;

void SteinerTreeProblem::dijkstra(vector<int> &u, vector<long long> &dist, vector<int> &p, vector<bool> &tree, priority_queue<pair<long long, int>, vector<pair<long long, int> >, greater<pair<long long, int> > > &pq, vector<bool> &I){
	if(p.empty()){
		long long INF = 1e15;
		p = vector<int>(n);
		dist = vector<long long>(n, INF);
		tree = vector<bool>(n, false);
	}

	for (int i = 0; i < u.size(); i++){
		tree[u[i]] = true;
		if(order.count(make_pair(dist[u[i]], u[i])) == 1) order.erase(make_pair(dist[u[i]], u[i])); 
		// Actualizamos order
		dist[u[i]] = 0; pq.push(make_pair(dist[u[i]], u[i]));
	}

	bool found = false;
	long long value = -1e15;
	while(!pq.empty()){
		pair<long long, int> tmp = pq.top(); pq.pop();
		long long d = tmp.first;
		int n = tmp.second;
		if(dist[n] != d) continue;
		FOREACH(vw, adj[n]){
			counter++;
			int v = vw->first;
			long long w = vw->second;
			if(!tree[v]){
				if(dist[v] > dist[n] + w){
					// Actualizamos order si se encontro algo mejor de algunos de los que aun no se han alcanzado
					if(I[v] && !tree[v]){
						if(order.count(make_pair(dist[v], v)) == 1) order.erase(make_pair(dist[v], v));
						order.insert(make_pair(dist[n] + w, v));
					}
					dist[v] = dist[n] + w;
					p[v] = n;
					pq.push(make_pair(dist[v], v));
					hsh_r[v] = 1;
				}
				else if(dist[v] == dist[n] + w && w != 0){
					if(rand()/(RAND_MAX + 1.0) < 1.0/hsh_r[v]) p[v] = n;
					hsh_r[v]++;
				}
			}
		}
		if (I[n] && (!tree[n])){
			//return;
			if(!found) value = dist[n], found = true;
		}
		if(found && dist[n] > value + currentTol) return;
	}
}

long long fitnessBest = 1e15;

/*
void SteinerTree::evaluateMinDistances(){
	long long best = 1e15;
	vector<bool> Ibest;
	int cnt_break = 0;
	vector<bool> bestTree;
	int threshold;
	if (random() % 2 == 0){
		threshold = 0;
	} else {
		threshold = random() % 5;
	}
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
			if (threshold != 0)
				threshold = 1 + random() % 4;
			long long mn = 1e15;
			u = -1;
			vector<int> ties;
			for(int i = 0; i < (SteinerTreeproblem->n); i++){
				if(I[i] && !tree[i] && dist[i] < mn){
					mn = dist[i];
				}
			}
			for(int i = 0; i < (SteinerTreeproblem->n); i++){
				if(I[i] && !tree[i] && dist[i] <= mn + threshold){
					ties.push_back(i);
				}
			}
			if (ties.size() == 0) break;
			u = ties[random() % (ties.size())];
			fitness += dist[u];
			vector<int> path;
			while(!tree[u]){
				path.push_back(u);
				cnt[u]++, cnt[p[u]]++;
				u = p[u];
			}
			reverse(path.begin(), path.end());
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
				sigset_t sign;
				sigemptyset(&sign);
				sigaddset(&sign, SIGTERM);
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
}*/


void SteinerTree::evaluateMinDistances(){
	long long best = 1e15;
	hsh_r = vector<int>(SteinerTreeproblem->n, 1);
	vector<bool> Ibest, ICompletebest;
	int cnt_break = 0;
	/*if(rand()%2 == 0) tol = 0;
	else tol = rand()%9 + 1;*/
	//int maxTol = generation / 100;
	//maxTol = min(maxTol, 5);
	currentTol = rand() % 2;
	if (currentTol != 0){
		currentTol = 1 + rand() % (maxThreshold + 1);
	}

	//int typeTransform = random() % 2;
	int selectionType = random() % 2;

	while(true){
		order.clear();
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
		int aux = 0;
		int connected = 1;

		/*
		vector<int> active;
		for(int i = 0; i < (SteinerTreeproblem->n); i++){
			if (I[i]){
				active.push_back(i);
			}
		}*/
		//cout << "Numero de activos: " << active.size() << endl;
		while(true){
     if (connected == SteinerTreeproblem->fs.size()){
        for (int i = 0; i < SteinerTreeproblem->n; i++){
          if ((I[i]) && (!tree[i])){
            if (SteinerTreeproblem->fixed[i]){
              cerr << "Error interno" << endl;
            }
            I[i] = false;
          }
        }
				break;
      }   

			aux++;
			long long mn = 1e15;
			u = -1;

			int cnt_r = 1;
			mn = order.begin()->first;
			long long t;
			FOREACH(v, order){
				if(v->first > mn + currentTol) break;
				if(rand()/(RAND_MAX + 1.0) < 1.0/cnt_r) u = v->second, t = v->first;
				cnt_r++;
			}

			if(u == -1) break;
			if (selectionType == 1){//Greedy
				u = -1;
				cnt_r = 1;
				FOREACH(v, order){
        	if(v->first > mn + currentTol) break;
        	if ((u != -1) && (!SteinerTreeproblem->fixed[u]) && (SteinerTreeproblem->fixed[v->second])){
          	cnt_r = 1;
          	u = v->second;
          	continue;
        	}
        	if ((u != -1) && (SteinerTreeproblem->fixed[u]) && (!SteinerTreeproblem->fixed[v->second])) continue;
        	if(rand()/(RAND_MAX + 1.0) < 1.0/cnt_r) u = v->second, t = v->first;
        	cnt_r++;
      	}

				if (!SteinerTreeproblem->fixed[u]){
 	       vector<int> nodes;
 	       FOREACH(v, order){
 	         if(v->first > mn + currentTol) break;
 	         nodes.push_back(v->second);
 	       }
 	       random_shuffle(nodes.begin(), nodes.end());
 	       vector < pair < int, pair<int, int > > > scoreAndNode;
 	       for (int i = 0; i < min((int)nodes.size(), 50); i++){
 	         int node = nodes[i];
 	         int score = 0;
 	         for (int i = 0; i < (SteinerTreeproblem->adj)[node].size(); i++){
 	           int otherNode = (SteinerTreeproblem->adj)[node][i].first;
 	           if ((SteinerTreeproblem->fixed[otherNode]) && ((!tree[otherNode]))){
 	             score++;
 	           }
 	         }
 	         scoreAndNode.push_back(make_pair(score, make_pair(random() % 100000, node)));
 	       }
 	       sort(scoreAndNode.begin(), scoreAndNode.end());
 	       u = scoreAndNode.back().second.second;
 	 	    }
			}

			if ((u != -1) && (SteinerTreeproblem->fixed[u])){
        connected++;
      }   

			order.erase(make_pair(dist[u], u)); // Eliminamos el ya utilizado
			//active[selectedIndex] = active.back();
			//active.pop_back();
			fitness += dist[u];
			vector<int> path;
			while(!tree[u]){
				path.push_back(u);
				cnt[u]++, cnt[p[u]]++;
				u = p[u];
				//printf("size = %d, u = %d, p = %d, it = %d\n", (int)path.size(), u, p[u], aux);
			}
			reverse(path.begin(), path.end());
			//for(int i = 0; i < (int)path.size(); i++) SteinerTreeproblem->dijkstra(path[i], dist, p, tree);
			//if (aux % 1000 == 0) cout << "Llama " << aux << endl;
			SteinerTreeproblem->dijkstra(path, dist, p, tree, pq, I);
			if (finished){
				printBest();
				exit(0);
			}
		}
		//cout << "El contador es: " << counter << endl;
		vector<bool> copyI = I;
		vector<bool> otherI(copyI.size());
		for(int u = 0; u < (SteinerTreeproblem->n); u++){
			if((SteinerTreeproblem->fixed)[u] || cnt[u] >= 1){
				I[u] = true;
			} else {
				I[u] = false;
			}
		}
		reset();
		I = copyI;

		if(best > fitness){
			//printf("best = %lld\n", best);                              
			best = fitness;
			Ibest = I;
			ICompletebest = IComplete;
			cnt_break = 0;
			if(fitness < Globalbest){
				cerr << "Error interno" << endl;
			}
		}
		else cnt_break++;
		if (cnt_break > 3) break; 
		int ones = 0;
		for(int u = 0; u < (SteinerTreeproblem->n); u++){
			if ((!(SteinerTreeproblem->fixed)[u]) && (I[u])){
				ones++;
			}
		}
		//cout<< "Hay " << ones << endl;
		int equals = 0;
		//if (typeTransform == 0){
			for(int u = 0; u < (SteinerTreeproblem->n); u++){
				if((SteinerTreeproblem->fixed)[u] || cnt[u] >= 3){
					if ((!(SteinerTreeproblem->fixed)[u]) && (I[u])) equals++;
					I[u] = true; 
				} else {
					I[u] = false;
				}
			}
		//}
	}
	I = Ibest;
	IComplete = ICompletebest;
	fitness = best;
	//printf("best = %lld, fitness = %lld\n", Globalbest, fitness);
	//printBest();
}

#ifndef __MA_H__
#define __MA_H__

struct ExtendedIndividual {
	SteinerTree ind;
	int dist;
};
extern int generation;
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
#include <sys/time.h>
#include <iostream>
#include <signal.h>

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
		if ((population.size() >= 8) && (population.size() % 2 == 0)){
			struct timeval currentTime;
			gettimeofday(&currentTime, NULL);
			double time2 = ((double) (currentTime.tv_sec) * 1.0e6 + (double) (currentTime.tv_usec))/1.0e6;
			double elapsed = time2 - time;
			if (elapsed > 30 * 60.0 / 500){
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
	int lastSelected;
	for (int i = parents.size(); i < N; i++){
		int first = getRandomInteger0_N(N - 1);
		int second = getRandomInteger0_N(N - 1);
		if (population[first]->ind.fitness <= population[second]->ind.fitness){
			/*if ((i % 2 == 1) && (first == lastSelected)){ 
				i--;
				continue;
			}*/
			lastSelected = first;
			parents.push_back(population[first]);
		} else {
			/*if ((i % 2 == 1) && (second == lastSelected)){
				i--;
				continue;
			}*/
			lastSelected = second;
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
	//if (generation > 100){
		for (int i = 0; i < offspring.size(); i++){
			offspring[i]->ind.dependentMutation(pm);
		}
	//}
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
	map<long long, int> acceptedFitness;
	acceptedFitness[population[0]->ind.fitness]++;
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
			//if ((acceptedFitness[all[indexBest]->ind.fitness] > N / 5.0) && (acceptedFitness[all[i]->ind.fitness] <= N / 5.0)){
			//	indexBest = i;
			//} else if ((acceptedFitness[all[indexBest]->ind.fitness] <= N/5.0) && (acceptedFitness[all[i]->ind.fitness] > N / 5.0)){
			//} else 
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
		acceptedFitness[all[indexBest]->ind.fitness]++;
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
	meanDistance /= ((population.size() * (population.size() - 1)) / 2);
	DI = meanDistance * 0.25;
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

#include <unistd.h>

int main(int argc, char **argv){
	//alarm(30 * 60 + 10);
	int N = 100;
	double pc = 1;
	double pm = 0.01;
	double finalTime = 30 * 60;
	MA ma(N, pc, pm, finalTime);
	//srand(time(NULL));
	srand(889);
	//string file = string(argv[1]);
	SteinerTreeProblem STP;
	SteinerTree::SteinerTreeproblem = &STP;
	ma.run();
}
