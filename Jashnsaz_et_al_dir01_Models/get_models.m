function Model = get_models()
% regulations are modeled as: 
% positive regulation (+1)
% negeative regulation (-1)
% none (0)
% get_ODE(model,input_profile) takes a given model and an input_profile,
% builds the ODE for that model and adds it to the model. 

% Model 1
model_index = 1;
Model{model_index}.n_nodes=4; % number of nodes of the model
Model{model_index}.A=[0, 0, 0, 0; % internode regulations matrix; 
                     0, 0, 0, 0; 
                    -1, 1, 0, 0; 
                     0, 0, 1, 0];

Model{model_index}.B=[-1,1,0,0]'; % regulations of model nodes via input 
Model{model_index}.C=[-1,-1,-1,-1]'; % basal regulations on the model nodes 
Model{model_index}.n_params = 2*nnz(Model{model_index}.A)+2*nnz(Model{model_index}.B)+2*nnz(Model{model_index}.C);

% Model 2
model_index = model_index + 1;
Model{model_index}.n_nodes=4;
Model{model_index}.A=[0, 0, 1, 0; 
                     0, 0, 0, 0; 
                    -1, 1, 0, 0; 
                     0, 0, 1, 0];

Model{model_index}.B=[-1,1,0,0]';
Model{model_index}.C=[-1,-1,-1,-1]';
Model{model_index}.n_params = 2*nnz(Model{model_index}.A)+2*nnz(Model{model_index}.B)+2*nnz(Model{model_index}.C);

% Model 3
model_index = model_index + 1;
Model{model_index}.n_nodes=4;
Model{model_index}.A=[0, 0, 1, 0; 
                     0, 0, 0,-1; 
                    -1, 1, 0, 0; 
                     0, 0, 1, 0];

Model{model_index}.B=[-1,1,0,0]';
Model{model_index}.C=[-1,-1,-1,-1]';
Model{model_index}.n_params = 2*nnz(Model{model_index}.A)+2*nnz(Model{model_index}.B)+2*nnz(Model{model_index}.C);

% Model 4
model_index = model_index + 1;
Model{model_index}.n_nodes=4;
Model{model_index}.A=[0, 0, 1, 1; 
                     0, 0, 0,-1; 
                    -1, 1, 0,-1; 
                     0, 0, 1, 0];

Model{model_index}.B=[-1,1,0,0]';
Model{model_index}.C=[-1,-1,-1,-1]';
Model{model_index}.n_params = 2*nnz(Model{model_index}.A)+2*nnz(Model{model_index}.B)+2*nnz(Model{model_index}.C);

%Model 5
model_index=model_index+1;
Model{model_index}.n_nodes=5;
Model{model_index}.A=[0, 0, 0, 1, 0; 
                     0, 0, 0, 0, 0; 
                     0, 0, 0, 0,-1; 
                    -1, 1, 1, 0, 0; 
                     0, 0, 0, 1, 0];
             
Model{model_index}.B=[-1,1,1,0,0]';
Model{model_index}.C=[-1,-1,-1,-1,-1]';
Model{model_index}.n_params = 2*nnz(Model{model_index}.A)+2*nnz(Model{model_index}.B)+2*nnz(Model{model_index}.C);

save('Models', 'Model');
end




