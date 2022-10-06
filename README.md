This repository contains the code for all algorithms discussed in the paper "Arcflow formulations and constraint generation frameworks for the two-bar charts packing problem" by Mathijs Barkel and Maxence Delorme. 

All algorithms are coded in C++. Moreover, note that all of our methods, except GA_LO, require the commercial solver Gurobi (we used version 9.5.0). The code is divided over 11 folders, each containing the code of one method. The different folders correspond to the following methods as named in Section 5 of our paper:
- 1_NAIVE       | ERZIN
- 2_CSP         | CSP
- 4_RELAXV1     | EUL-REL
- 5_RELAXV2     | LINK-REL
- 6_VPP         | 2DV
- 8_EXACTV1     | EUL-FLOW
- A_EXACTV2     | LINK-FLOW
- B_LP          | continous relaxation of ERZIN
- C_LPV1        | continuous relaxation of EUL-REL
- D_LPV2        | continous relaxation of LINK-REL
- E_GALO        | GA_LO

Each folder contains the same substructure. For example, 8_EXACTV1 contains the following subfolders:
- EXACTV1.cpp: the main back-end code for the method EUL-FLOW
- EXACTV1.h: the header file corresponding to EXACTV1.cpp
- helper_functions.cpp: contains a number of secondary functions (this file is the same for each of the 11 main folders)
- helper_functions.h: the header file corresponding to helper_functions.cpp (this file is the same for each of the 11 main folders)
- main.cpp: the front-end code for using the method EUL-FLOW
- main.h: the header file corresponding to EXACTV1.cpp
- makefile: used for compiling under linux (it needs to be updated by the user)

Once compiled, the following command can be used to run the algorithm:
  ./PROGRAM "./PATH_INSTANCE" "NAME_INSTANCE" "./PATH_AND_NAME_OUTPUT_GENERAL" 
where
- PROGRAM is the name of the compiled software (e.g., EXACTV1)
- ./PATH_INSTANCE is the relative path of the folder where the instance to solve is located
- NAME_INSTANCE is the name of the instance to solve
- ./PATH_AND_NAME_OUTPUT_GENERAL is the name of the file (together with its relative path) where performance metrics (such as the optimality status, the CPU time required, or the number of variables) are stored after solving an instance

Moreover, "Instances.rar" contains a txt-file for each of our test instances. There are 6 main folders, each corresponding to a different instance type (see Section 5 of the paper):
- 1_U-GEN 	  (C1)
- 2_U-SMA	    (CS)
- 3_U-MED	    (CM)
- 4_U-BIG	    (CB)
- 5_TRIPLETS	(CT)
- 6_DONUTS	  (CD)

The name of each txt-file is of the following format: "A_B_C_D_E(_F_G)_H", where:
- A: a code referring to the instance type (C1/CS/CM/CB/CT/CD)
- B: the total number of bar charts (Sigma) 
- C: the bin capacity (c)
- D and E: refer to the way that item weights are generated. For example, all U-SMA instances have D = 0 and E = 0.1, meaning that each item weight is uniformly distributed in [1, 0.1c]. 
- F (only for DONUTS instances): the number of donuts (either 1 or 3)
- G (only for DONUTS instances): refers to the number of bins per donut (z2-bar) (either 2 or 0.2, where 0.2 means 0.2 * z1-bar)
- H: the instance index (there are 10 instances for each combination of parameters)

Each txt-file is structured as follows:
- the first line contains the number of distinct bar chart types (n)
- the second line contains the bin capacity (c)
- the remaining (n) lines all contain a triplet of numbers. Each triplet refers to a unique bar chart type, giving respectively:
    - the weight of the first item (w_i1)
    - the weight of the second item (w_i2)
    - the number of bar charts of that type (q_i)
