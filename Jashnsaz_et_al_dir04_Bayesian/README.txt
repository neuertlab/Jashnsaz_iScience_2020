These codes implement Bayesian analysis. (Figures 4J-4L).

1. run_Bayes.m specifies local computer versus cluster, kinetic types (training datasets IDs), model free parameters, and the number of chains for each condition. It runs the function Bayesian.m for each condition. 

2. Bayesian.m imports the true model, the true parameters, and the simulated data from the directory called Models_TrueParams_simData. The simulated data (sim_data.mat) contains simulated data for all 58 input profiles given in Figure 1H. 
It gets training data according to the IDs specified as train_data_IDs in run_Bayes.m. 

3. Bayesian.m uses the following functions: 
Get_ODE.m to construct model ODE for all train and test data input profiles. 
get_simp_OBJ.m to build the objective function. 
get_predictions.m and get_simp_sol.m to calculate fits (training datasets) and predictions (remaining datasets) for input profiles from all the 58 datasets for 100 parameters sets drawn from the posterior.  

4. jobs_manager_local.sh is used to submit jobs on local computer. 
after rooting to the directory, using commands: 
chmod u+x jobs_manager_local.sh
./jobs_manager_local.sh
which executes run_jobs.m, within which the jobs could be specified. 

4. jobs_manager_cluster.slurm is used to submit jobs on a cluster using slurm. 
It submits all the 30 jobs specified in run_Bayes.m. 
 

