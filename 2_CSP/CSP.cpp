#include "CSP.h"

Solution CSP(const Instance& inst) {

	// local variable
	double start = getCPUTime(); 			// starting time
	Solution sol; sol.Ncuts = 0;				// solution to return
	
	// Transform instance
	vector<int> items (inst.c + 1,0);
	for (int i = 0; i < inst.m; i++){
		items[inst.BCs[i][0]]+=inst.BCs[i][2];
		items[inst.BCs[i][1]]+=inst.BCs[i][2];
	}
	
	// Create graph
	vector<bool> hasBeenVisited (inst.c + 1,false); hasBeenVisited[0] = true;
	vector<vector<int>> arcs;
	for (int i = inst.c ; i >= 0; i--){
		vector<bool> hasBeenVisitedHere(inst.c + 1,false);
		for (int j = 0 ; j < items[i]; j++){
			for (int k = inst.c - i ; k >= 0; k--){
				if(hasBeenVisited[k] && !hasBeenVisitedHere[k]){
					hasBeenVisitedHere[k] = true;
					arcs.push_back({k,k+i,i});
					hasBeenVisited[k+i] = true;
				}
			}
		}
	}
	
	// Add loss arcs
	vector<int> visitiedNodes; 
	for(int i =0;i<inst.c;i++) if (hasBeenVisited[i]) visitiedNodes.push_back(i);
	visitiedNodes.push_back(inst.c);
	for(int i =0;i<visitiedNodes.size()-1;i++) arcs.push_back({visitiedNodes[i],visitiedNodes[i+1],-1});
		
    // create a model
    GRBEnv env = GRBEnv();              	// create an environment
    GRBModel model = GRBModel(env);         // create a new model

    // declaration of the variables for the model
    vector<GRBVar> fa(arcs.size());
	
    // initizalization of the variables for the model
    for (int i = 0; i < arcs.size(); i++) { 
		cout << arcs[i][0] << " " << arcs[i][1] << " " << arcs[i][2] << endl;
		if(arcs[i][2] != -1) fa[i] = model.addVar(0, items[arcs[i][2]], 0, GRB_INTEGER);
		else fa[i] = model.addVar(0, inst.UB, 0, GRB_INTEGER);
    }
    model.update();

	// declare linear expressions
    vector<GRBLinExpr> fIn (inst.c+1, 0);    		// the amount of flow entering each vertex
    vector<GRBLinExpr> fOut (inst.c+1, 0);   		// the amount of flow leaving each vertex
    vector<GRBLinExpr> typeUsed(items.size(), 0); 	// the amount of arcs used of each BC type 
    
	// calculate the linear expressions	
    for (int i = 0; i < arcs.size(); i++) {   				// loop over all arcs
        fIn[arcs[i][1]] += fa[i];        					// inflow
        fOut[arcs[i][0]] += fa[i];       					// outflow
		if(arcs[i][2] >= 0) typeUsed[arcs[i][2]] += fa[i];	// number of items used of certain type
    }
	model.update();

    // create flow conservation constraints 
    for (int v = 1; v < inst.c; v++) {       	// loop over all vertices
		if (hasBeenVisited[v])
			model.addConstr(fIn[v] == fOut[v]); // inflow = outflow
    }
	model.addConstr(fOut[0] == fIn[inst.c]);


    // create item type constraints
    for (int i = 0; i < items.size(); i++) {    // loop over all item types
		if (items[i] >= 1){
			if (items[i] == 1) model.addConstr(typeUsed[i] == items[i]); // demand met exactly if 1
			else model.addConstr(typeUsed[i] >= items[i]); 		  		 // demand met otherwise
		}
    }

    // set the objective: minimize the number of flow going out of 0
    model.setObjective(fOut[0], GRB_MINIMIZE);

    // change some settings
    model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
    model.getEnv().set(GRB_IntParam_Threads, 1);
    model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods
    model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);

    // find the optimal solution
	sol.timeP = getCPUTime() - start;			// preprocessing time
    model.optimize();

    // store the results in a Solution object
    sol.Nvar = model.get(GRB_IntAttr_NumVars);       // number of variables
    sol.Nconstr = model.get(GRB_IntAttr_NumConstrs); // number of constraints
    sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);      // number of non-zero coefficients		
	sol.opt = 0;
	sol.LB = ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
		
    if (model.get(GRB_IntAttr_SolCount) >= 1) { 		// if a solution has been found
        sol.UB = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		
		// get bins
		int binIdx = 0;
		vector<vector<vector<int>>> arcsUsed (inst.c+1);
		int nbRem = 0;
		for (int i = 0; i < arcs.size(); i++) {
			for (int k = 0; k < ceil(fa[i].get(GRB_DoubleAttr_X) - EPSILON); k++){  
				arcsUsed[arcs[i][0]].push_back(arcs[i]);
				nbRem++;
			}
		}

		while (arcsUsed[0].size() > 0){ 
			vector<int> bin; int tail = 0; int head = 0; int load = 0;
			while(tail != inst.c){ 
				head = arcsUsed[tail].back()[1];
				if(arcsUsed[tail].back()[2] >= 0){
					bin.push_back(arcsUsed[tail].back()[2]);
					load += arcsUsed[tail].back()[2];
				}
				arcsUsed[tail].pop_back(); 
				tail = head;
				nbRem--;
			}
			cout << "Bin " << binIdx << "\t" << load << "\t ["; binIdx++;
			for (int j = 0; j <  bin.size(); j++) cout << bin[j] << " ";
			cout << "]" << endl;
		}
    }

    sol.timeT = getCPUTime() - start;	
    return sol;
}

