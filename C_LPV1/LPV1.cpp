#include "LPV1.h"

Graph makeGraph1(const Instance& inst) {

	// Local variables
	Graph G;
	vector<bool> Sbool(inst.c + 1, false); // set of reachable start vertices S
    
	// Initialisation  
	G.NTI.resize(inst.c + 1, vector<int>(inst.c + 1, -1));
    Sbool[0] = true;
 
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

    // Find all reachable vertices by 'shifting' the base layer
	int nbArcs = G.A.size();							// fix the number of arcs for all iteration
    for (int s = 1; s < inst.c+1; s++){
		if (Sbool[s]){
			for(int a =0; a<nbArcs;a++){
				if (G.ITN[G.A[a][1]][0] + s <= inst.c) {
					newC1 = G.ITN[G.A[a][0]][0] + s;	
					newC2 = G.ITN[G.A[a][0]][1];
					if (G.NTI[newC1][newC2] == -1) {	// if the new tail does not exist yet, create it
						G.NTI[newC1][newC2] = G.Nvert;	
						G.ITN.push_back({newC1,newC2});
						G.Nvert++; 		
					}
					newC3 = G.ITN[G.A[a][1]][0] + s;	
					newC4 = G.ITN[G.A[a][1]][1];
					if (G.NTI[newC3][newC4] == -1) {	// if the new head does not exist yet, create it
						G.NTI[newC3][newC4] = G.Nvert;	
						G.ITN.push_back({newC3,newC4});
						G.Nvert++; 		
					}
					G.A.push_back({G.NTI[newC1][newC2],G.NTI[newC3][newC4],G.A[a][2]});   
				}
			}
		}
    }

    // Add the loss arcs 
    for (int i = 1; i < G.Nvert; i++) {     	
		newC1 = G.ITN[i][1];	
		newC2 = 0;
		if (G.NTI[newC1][newC2] == -1) {	// if the new head does not exist yet, create it
			G.NTI[newC1][newC2] = G.Nvert;	
			G.ITN.push_back({newC1,newC2});
			G.Nvert++; 		
		}
		G.A.push_back({i, G.NTI[newC1][newC2], -1}); 
    }
	
	// Erase duplicates in arcs
	set<vector<int>> s(G.A.begin(),G.A.end());
	G.A.assign(s.begin(),s.end());
	G.Narcs = G.A.size();
	cout << "---------------------------------------" << endl;
	
    // return the Graph object
    return G;
}

Solution LPV1(const Instance& inst) {

	// local variable
	double start = getCPUTime(); 			// starting time
	Solution sol; 							// solution to return
	Graph G = makeGraph1(inst);   
	
    // create a model
    GRBEnv env = GRBEnv();              	// create an environment
    GRBModel model = GRBModel(env);         // create a new model

    // declaration of the variables for the model
    vector<GRBVar> fa(G.Narcs);
	
    // initizalization of the variables for the model
    for (int i = 0; i < G.Narcs; i++) { 
		// cout << G.ITN[G.A[i][0]][0] << "," << G.ITN[G.A[i][0]][1] << " " << G.ITN[G.A[i][1]][0] << "," << G.ITN[G.A[i][1]][1] << " " << G.A[i][2] << endl;
		if(G.A[i][2] != -1) fa[i] = model.addVar(0, inst.BCs[G.A[i][2]][2], 0, GRB_CONTINUOUS);
		else fa[i] = model.addVar(0, inst.UB, 0, GRB_CONTINUOUS);
    }
    model.update();

	// declare linear expressions
    vector<GRBLinExpr> fIn (G.Nvert, 0);    		// the amount of flow entering each vertex
    vector<GRBLinExpr> fOut (G.Nvert, 0);   		// the amount of flow leaving each vertex
    vector<GRBLinExpr> typeUsed(inst.m + 1, 0); 	// the amount of arcs used of each BC type 
    
	// calculate the linear expressions	
    for (int i = 0; i < G.Narcs; i++) {   		// loop over all arcs
        fIn[G.A[i][1]] += fa[i];        		// inflow
        fOut[G.A[i][0]] += fa[i];       		// outflow
        typeUsed[G.A[i][2]+1] += fa[i];      	// number of BCs used of certain type
    }
	model.update();

    // create flow conservation constraints 
    for (int v = 0; v < G.Nvert; v++) {       	// loop over all vertices
        model.addConstr(fIn[v] == fOut[v]); 	// inflow = outflow
    }

    // create constraint ensuring at least one cycle containing (0, 0)
    model.addConstr(fOut[0] >= 1); 

    // create BC type constraints
    for (int i = 0; i < inst.m; i++) {           // loop over all BC types
		if (inst.BCs[i][2] == 1) model.addConstr(typeUsed[i + 1] == inst.BCs[i][2]); 	// demand met exactly if 1
		else model.addConstr(typeUsed[i + 1] >= inst.BCs[i][2]); 						// demand met otherwise
    }

    // set the objective: minimize the number of loss arcs used
    model.setObjective(typeUsed[0], GRB_MINIMIZE);

    // change some settings
    model.getEnv().set(GRB_IntParam_Threads, 1);
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

