#include "VPP.h"

Graph makeGraphVP(const Instance& inst) {

	// Local variables
	Graph G;
	vector<bool> Sbool(inst.c + 1, false); // set of reachable start vertices S
    
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
						Sbool[newC4] = true;
					}
					G.A.push_back({G.NTI[newC1][newC2],G.NTI[newC3][newC4],k});   
				}
				else break;
			}
		}
    }
	
	// Erase duplicates in arcs
	set<vector<int>> s(G.A.begin(),G.A.end());
	G.A.assign(s.begin(),s.end());
	G.Narcs = G.A.size();
	cout << "---------------------------------------" << endl;
	
    // return the Graph object
    return G;
}

Solution VPP(const Instance& inst) {

	// local variable
	double start = getCPUTime(); 			// starting time
	Solution sol; sol.Ncuts = 0;			// solution to return
	Graph G = makeGraphVP(inst);   
	
    // create a model
    GRBEnv env = GRBEnv();              	// create an environment
    GRBModel model = GRBModel(env);         // create a new model

    // declaration of the variables for the model
    vector<GRBVar> fa(G.Narcs);
	
    // initizalization of the variables for the model
    for (int i = 0; i < G.Narcs; i++) { 
		// cout << G.ITN[G.A[i][0]][0] << "," << G.ITN[G.A[i][0]][1] << " " << G.ITN[G.A[i][1]][0] << "," << G.ITN[G.A[i][1]][1] << " " << G.A[i][2] << endl;
		fa[i] = model.addVar(0, inst.BCs[G.A[i][2]][2], 0, GRB_INTEGER);
    }
    model.update();

	// declare linear expressions
    vector<GRBLinExpr> fIn (G.Nvert, 0);    	// the amount of flow entering each vertex
    vector<GRBLinExpr> fOut (G.Nvert, 0);   	// the amount of flow leaving each vertex
    vector<GRBLinExpr> typeUsed(inst.m, 0); 	// the amount of arcs used of each BC type 
    
	// calculate the linear expressions	
    for (int i = 0; i < G.Narcs; i++) {   		// loop over all arcs
        fIn[G.A[i][1]] += fa[i];        		// inflow
        fOut[G.A[i][0]] += fa[i];       		// outflow
        typeUsed[G.A[i][2]] += fa[i];      		// number of BCs used of certain type
    }
	model.update();

    // create flow conservation constraints 
    for (int v = 1; v < G.Nvert; v++) {       	// loop over all vertices
        model.addConstr(fIn[v] >= fOut[v]); 	// inflow >= outflow
    }

    // create BC type constraints
    for (int i = 0; i < inst.m; i++) {           // loop over all BC types
		if (inst.BCs[i][2] == 1) model.addConstr(typeUsed[i] == inst.BCs[i][2]); 	// demand met exactly if 1
		else model.addConstr(typeUsed[i] >= inst.BCs[i][2]); 						// demand met otherwise
    }

    // set the objective: minimize the number of loss arcs used
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
	sol.LB = 2*ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
		
    if (model.get(GRB_IntAttr_SolCount) >= 1) { 		// if a solution has been found
        sol.UB = 2*ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		
		// get bin for each bar chart
		vector<vector<vector<int>>> arcsUsed (G.Nvert);
		vector<vector<int>> bins;
		for (int i = 0; i < G.Narcs; i++) {
			for (int k = 0; k < ceil(fa[i].get(GRB_DoubleAttr_X) - EPSILON); k++){  
				arcsUsed[G.A[i][0]].push_back(G.A[i]);
			}
		}

		while (arcsUsed[0].size() > 0){
			vector<int> bin;
			int tail = 0; int head; 
			// Create bin
			while(arcsUsed[tail].size() > 0){
				head = arcsUsed[tail].back()[1];
				bin.push_back(arcsUsed[tail].back()[2]); 
				arcsUsed[tail].pop_back(); 
				tail = head;
			}
			bins.push_back(bin);
		}

		// get bin for each BC type
		sol.assignedBin.resize(inst.m);
		for(int i = 0; i < bins.size();i++){
			for(int j = 0; j < bins[i].size();j++){
				sol.assignedBin[bins[i][j]].push_back(2*i);
			}
		}
    }

    sol.timeT = getCPUTime() - start;	
    return sol;
}

