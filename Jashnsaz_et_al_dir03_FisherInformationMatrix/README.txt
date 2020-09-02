These codes implement Fisher Information Matrix analysis. 

1. optimal_experiment_design.m imports the true model, the true parameters, and the simulated data from the directory called Models_TrueParams_simData. 

2. optimal_experiment_design.m calculates individual FIM matrixes (or imports them FIMs.mat) 
for datasets 1 to 36 given in Figure 1H using the functions get_FIM.m, Get_ODE.m, and Get_Sens_ODE.m. 

3. optimal_experiment_design.m then considers all combinatorial experiments containing n (for n=1,2,...6) datasets out of 36 datasets. For each condition, it calculates and outputs the experiment IDs for the D-Optimal experiment  design, that is the datasets minimizing the determinant of inverse FIM (FIM^{-1}). 

4. optimal_experiment_design.m plots determinant(FIM^{-1}) versus n for D-Optimal experiment design compared to STEPS and DIVERSE kinetics. Figure 4H. 


