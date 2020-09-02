function run_Bayes(nodeID)

local_vs_cluster_options  = {'local','cluster'}; 
local_vs_cluser = local_vs_cluster_options{1}; % to specify machine type

kinetics_names = {'STEPS','DIVERSE07M', 'DOPTIMAL'}; % kinetics types
DATA_IDs = [[1:6]; [6:6:36]; [5 6 12 23 35 36]]; % Refer to Figure 1H in the manuscript

free_model_parameters = [1:22]; % full model parameters are free
n_chains = 10; % number of independent chains for each conditions

i=1; 
for k=1:length(kinetics_names)
    for chain=1:n_chains

        variables{i}.kinetics_name=kinetics_names{k};
        variables{i}.train_data_IDs = DATA_IDs(k,:); 
        variables{i}.chain = chain;
        i = i + 1; 

    end
end

Bayesian(local_vs_cluser,free_model_parameters,variables{nodeID}.kinetics_name,variables{nodeID}.train_data_IDs,variables{nodeID}.chain); 
end