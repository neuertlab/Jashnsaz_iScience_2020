These codes generate the models and visualize them. 

1. get_models.m generates the 5 models used in Figure S5. Model3 is the true model. 
get_models.m gets the regulations matrixes (Figure S1B) from the user. Saves Models.mat 

2. visualize_models.m imports the Models.mat and uses the customized network_plot.m 
function to visualize all the models network diagrams (similar to Figure S1A). 
It outputs Models.png for all (5) models and ModelsGraphs directory containing the graphs
for the individual models.


