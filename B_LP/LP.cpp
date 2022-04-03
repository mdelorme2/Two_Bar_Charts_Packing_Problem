#include "LP.h"

Solution LP(const Instance& inst) {

	// local variable
	int n = inst.UB; 						// upper bound on the number of necessary bins
	double start = getCPUTime(); 			// starting time
    Solution sol; 							// solution to return
	
    // create a model
    GRBEnv env = GRBEnv();              	// create an environment
    GRBModel model = GRBModel(env);         // create a new model

    // declaration of the variables for the model
    vector<vector<GRBVar> > x;				// xij = number of times BC type i is placed in cell j
    x.resize(inst.m, vector<GRBVar>(n)); 
    vector<GRBVar> y(n); 					// yj = 1 if cell j is used, 0 otherwise

    // initizalization of the variables for the model
    for (int j = 0; j < n; j++) {       	
        y[j] = model.addVar(0, 1, 0, GRB_CONTINUOUS);
        for (int i = 0; i < inst.m; i++) {   	
            x[i][j] = model.addVar(0, inst.BCs[i][2], 0, GRB_CONTINUOUS);
        }
    }

    model.update();

    // create linear expressions
    GRBLinExpr obj = 0;                 	// the total number of cells used
    vector<GRBLinExpr> assigned(inst.m, 0); // the amount of cells that a BC type is assigned to
    vector<GRBLinExpr> height(n, 0);   		// the sum of all bar heights placed in each cell
    for (int j = 0; j < n - 1; j++) { 		// for all cells excepted the last one
        obj += y[j]; 						// the objective is the sum of all y-variables
        for (int i = 0; i < inst.m; i++) { 		// for all BC types
            height[j] 	+= inst.BCs[i][0] * x[i][j];
			height[j+1] += inst.BCs[i][1] * x[i][j];
            assigned[i] += x[i][j];
        }
    }
    obj += y[n - 1];						// add the last cell
    
    model.update();

    // create assignment constraints
    for (int i = 0; i < inst.m; i++) {      										// loop over all BC types
		if (inst.BCs[i][2] == 1) model.addConstr(assigned[i] == inst.BCs[i][2]); 	// this BC type must meet the demand exactly if the demand is 1
		else model.addConstr(assigned[i] >= inst.BCs[i][2]); 						// this BC type must just meet the demand otherwise
    }

    // create cell height constraints and symmetry breaking constraints
    for (int j = 0; j < n-1; j++) {           		// loop over all cells (except the last)
        model.addConstr(height[j] <= inst.c*y[j]); 	// the total height shouldn't exceed the cell's capacity
        model.addConstr(y[j] >= y[j + 1]);    		// the next cell can only be used if this cell is used too
    }
    model.addConstr(height[n-1] <= inst.c * y[n-1]); // special case for last cell

    // set the objective: minimize obj (= sum of all y-variables)
    model.setObjective(obj, GRB_MINIMIZE);

    // change some settings
    model.getEnv().set(GRB_IntParam_Threads, 1);
    model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods

    // find the optimal solution
	sol.timeP = getCPUTime() - start;			// preprocessing time
	
    // store the results in a Solution object
    if(inst.m * n >5000000){
		model.getEnv().set(GRB_DoubleParam_TimeLimit, 5);
	}
	model.optimize();
    sol.LP = model.get(GRB_DoubleAttr_ObjVal);
	sol.Nvar = model.get(GRB_IntAttr_NumVars);       // number of variables
    sol.Nconstr = model.get(GRB_IntAttr_NumConstrs); // number of constraints
    sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);      // number of non-zero coefficients		

    sol.timeT = getCPUTime() - start;	
    return sol;
}

