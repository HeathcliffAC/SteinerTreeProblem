class SteinerTreeProblem{
	public:
		SteinerTreeProblem(const string &fileName);
		SteinerTreeProblem(){
		
		}
		~SteinerTreeProblem(){
		
		}
		vector<vector<int> > adj;
		vector<pair<int, int> > edges;
};

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
		vector<bool> I;
		long long fitness;
		SteinerTreeProblem *SteinerTreeproblem;
};
