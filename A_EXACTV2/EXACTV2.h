#ifndef NAIVE_H
#define NAIVE_H

using namespace std;
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <math.h> 
#include "gurobi_c++.h"
#include "helper_functions.h"

struct Graph
{
    vector<vector<int>> A;
    int Narcs, Nvert;
	vector<vector<int>> ITN;
	vector<vector<int>> NTI;
};

class mycallbackT2: public GRBCallback
{
	public:
		int Ncuts, UB;
		Graph G;
		vector<GRBVar> fa;
		vector<GRBVar> la;
		vector<vector<int>> lossArcs;
	    
		mycallbackT2(const Graph& xG, const vector<GRBVar>& xfa, const vector<GRBVar>& xla, const vector<vector<int>>& xlossArcs, const int& xUB);
		int getNcuts();

	protected:
		void callback();
};

Solution EXACTV2(const Instance& inst);
Graph makeGraph2(const Instance& inst);

#endif 