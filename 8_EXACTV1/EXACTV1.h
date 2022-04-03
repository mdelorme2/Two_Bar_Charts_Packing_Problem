#ifndef EXACTV1C2_H
#define EXACTV1C2_H

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

class mycallbackT1: public GRBCallback
{
	public:
		int Ncuts;
		Graph G;
		vector<GRBVar> fa;
	    vector<int> demand;
		int UB;
		
		mycallbackT1(const Graph& xG, const vector<GRBVar>& xfa, const vector<int>& xdemand, const int& xUB);
		int getNcuts();

	protected:
		void callback();
};

Solution EXACTV1(const Instance& inst);
Graph makeGraph1(const Instance& inst);

#endif 