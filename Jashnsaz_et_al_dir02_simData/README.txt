These codes use the true model and the true parameters, and generate simulated data for 54 input profiles given in Figure 1H for pathway (HOG) activation and visualize them. 

1. get_simData.m imports the true model and the true parameters from the directory called Models_and_TrueParams. 
It generates input profiles over different kinetic types and intensities.  
It  uses Get_ODE.m to build the ODEs for the true model upon each input profiles, 
which are then used to generate the simulated data sets. 
It saves the data, and plots and saves the inputs and the responses into simData directory. 
(Such datasets are given in Figure S2). 


