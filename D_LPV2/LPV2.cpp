#include "LPV2.h"

Graph makeGraph2(const Instance& inst) {
	
	// Local variables
	Graph G;
    
	// Initialisation  
	G.NTI.resize(inst.c + 1, vector<int>(inst.c + 1, -1));
 
    int newC1, newC2, newC3, newC4;
	 
	// First node (0,0) is active
    G.NTI[0][0] = 0;               		
	G.ITN.push_back({0,0}); 
    G.Nvert = 1;                   		

	// Loop over all BCs
    for (int k = 0; k < inst.m; k++) {  	
 		int nbNodes = G.Nvert;              			// fix the number of nodes for that iteration
		vector<bool> alreadyVisited (nbNodes,false);
		for (int i = 0; i < nbNodes; i++) {     		// loop over all the existing vertices
			for (int j = 0; j < inst.BCs[k][2]; j++) { 	// for each unit of demand of BC type k
				newC1 = G.ITN[i][0] + j * inst.BCs[k][0];	
				newC2 = G.ITN[i][1] + j * inst.BCs[k][1];
				newC3 = G.ITN[i][0] + (j+1) * inst.BCs[k][0];	
				newC4 = G.ITN[i][1] + (j+1) * inst.BCs[k][1];
				int idxN = G.NTI[newC1][newC2];
				if(idxN >= 0 && idxN < nbNodes){
					if (alreadyVisited[idxN] == true) break;
					else alreadyVisited[idxN] = true;
				}
				if(newC3 <= inst.c && newC4 <= inst.c){ // if one unit of the BC type can be included
					if (G.NTI[newC3][newC4] == -1) {	// if the new head does not exist yet, create it
						G.NTI[newC3][newC4] = G.Nvert;	
						G.ITN.push_back({newC3,newC4});
						G.Nvert++; 			
					}
					G.A.push_back({G.NTI[newC1][newC2],G.NTI[newC3][newC4],k});   
				}
				else break;
			}
		}
    }
	
	// Add transition arcs
    for (int i = 0; i < G.Nvert; i++) G.A.push_back({i, 0, -1}); 
	
	// Erase duplicates in arcs
	set<vector<int>> s(G.A.begin(),G.A.end());
	G.A.assign(s.begin(),s.end());
	G.Narcs = G.A.size();
	cout << "---------------------------------------" << endl;
	
    // return the Graph object
    return G;
}

Solution LPV2(const Instance& inst) {

	// local variable
	double start = getCPUTime(); 				// starting time
	Solution sol; 								// solution to return
	Graph G = makeGraph2(inst);   
	vector<bool> activeNodes (inst.c+1,false);	// for loss arcs
	
    // create a model
    GRBEnv env = GRBEnv();              	// create an environment
    GRBModel model = GRBModel(env);         // create a new model

    // declaration of the variables for the model
    vector<GRBVar> fa(G.Narcs);
	
    // initizalization of the variables for the model
    for (int i = 0; i < G.Narcs; i++) { 
		activeNodes[inst.c - G.ITN[G.A[i][0]][0]] = true; activeNodes[G.ITN[G.A[i][0]][1]] = true;
	//	cout << G.ITN[G.A[i][0]][0] << "," << G.ITN[G.A[i][0]][1] << " " << G.ITN[G.A[i][1]][0] << "," << G.ITN[G.A[i][1]][1] << " " << G.A[i][2] << endl;
		if(G.A[i][2] != -1) fa[i] = model.addVar(0, inst.BCs[G.A[i][2]][2], 0, GRB_CONTINUOUS);
		else fa[i] = model.addVar(0, inst.UB, 0, GRB_CONTINUOUS);
    }
	
	// create loss arcs structure
	vector<vector<int>> lossArcs;
	activeNodes[inst.c] = true;
	vector<int> isActive; for(int i = 0; i < inst.c+1;i++) if (activeNodes[i]) isActive.push_back(i);
	for (int i = 0; i < isActive.size()-1;i++) lossArcs.push_back({isActive[i],isActive[i+1]});
	vector<GRBVar> la(lossArcs.size());
	for (int i = 0; i < lossArcs.size(); i++){
		// cout << lossArcs[i][0] << " " << lossArcs[i][1] << endl;
		la[i] = model.addVar(0, inst.UB, 0, GRB_CONTINUOUS);
	}
    model.update();	

	// declare linear expressions
    vector<GRBLinExpr> fIn (G.Nvert, 0);    		// the amount of flow entering each vertex
    vector<GRBLinExpr> fOut (G.Nvert, 0);   		// the amount of flow leaving each vertex
    vector<GRBLinExpr> tIn (inst.c+1, 0);    		// the amount of flow entering each transition vertex
    vector<GRBLinExpr> tOut (inst.c+1, 0);   		// the amount of flow leaving each transition vertex
    vector<GRBLinExpr> typeUsed(inst.m + 1, 0); 	// the amount of arcs used of each BC type 
    
	// calculate the linear expressions	
    for (int i = 0; i < G.Narcs; i++) {   		// loop over all arcs
        fIn[G.A[i][1]] += fa[i];        		// inflow
        fOut[G.A[i][0]] += fa[i];       		// outflow
		if (G.A[i][2] == -1){
			tOut[inst.c - G.ITN[G.A[i][0]][0]] += fa[i];
			tIn[G.ITN[G.A[i][0]][1]] += fa[i];
		}
        typeUsed[G.A[i][2]+1] += fa[i];    		// number of BCs used of certain type
    }
	for (int i = 0; i < lossArcs.size(); i++){
		tOut[lossArcs[i][0]] += la[i];
		tIn[lossArcs[i][1]] += la[i];
	}
	
	model.update();

    // create flow conservation constraints 
    for (int v = 0; v < G.Nvert; v++) {       	// loop over all vertices
        model.addConstr(fIn[v] == fOut[v]); 	// inflow = outflow
    }
    for (int v = 0; v < isActive.size(); v++) { // loop over all transition vertices	
		model.addConstr(tIn[isActive[v]] == tOut[isActive[v]]);
	}

    // create constraint ensuring at least one cycle containing (0, 0)
    model.addConstr(fOut[0] >= 1); model.addConstr(tOut[0] >= 1); 

    // create BC type constraints
    for (int i = 0; i < inst.m; i++) {           // loop over all BC types
		if (inst.BCs[i][2] == 1) model.addConstr(typeUsed[i + 1] == inst.BCs[i][2]); 	// demand met exactly if 1
		else model.addConstr(typeUsed[i + 1] >= inst.BCs[i][2]); 						// demand met otherwise
    }

    // set the objective: minimize the number of loss arcs used
    model.setObjective(typeUsed[0], GRB_MINIMIZE);

    // change some settings
    model.getEnv().set(GRB_IntParam_Threads, 1);
//	model.getEnv().set(GRB_IntParam_Crossover, 0); 
    model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods

    // find the optimal solution
	sol.timeP = getCPUTime() - start;			// preprocessing time
    model.optimize();

    // store the results in a Solution object
    sol.Nvar = model.get(GRB_IntAttr_NumVars);       // number of variables
    sol.Nconstr = model.get(GRB_IntAttr_NumConstrs); // number of constraints
    sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);      // number of non-zero coefficients		
    sol.LP = model.get(GRB_DoubleAttr_ObjVal);

    sol.timeT = getCPUTime() - start;	
    return sol;
}

