#ifndef LPV1_H
#define LPV1_H

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

Solution LPV1(const Instance& inst);
Graph makeGraph1(const Instance& inst);

#endif 