function model = Get_Sens_ODE(model,input_profile,parameter)
%%% Get the ODEs for a model, along with the sensitivities for the parameter
% specified by parameter. 
% get the number of nodes, make empty ODE. 

n_nodes = size(model.A,1);
tmp_ODE_A = sym(zeros(n_nodes,1));
tmp_ODE_B = sym(zeros(n_nodes,1));
% get the number of parameters and 
% specify parameter vector for the models. 
n_parameters = 2*nnz(model.A)+2*nnz(model.B)+2*nnz(model.C);
sym_K = sym('K',[n_parameters 1]);
sym_X = sym('X',[n_nodes 1] );
param_count = 1;
% symbolic time
t = sym('t');
% use a seperate symbol for basal deactivation terms. 
d = sym('d');
% symbolic input (u(t)) for symbolic ODE
s = sym('s');

for j=1:n_nodes
    % add the terms for inputs regulations
    if model.B(j) == 1
        tmp_ODE_B(j) = tmp_ODE_B(j) + ...
            sym_K(param_count)*(1-sym_X(j))/((1-sym_X(j))+sym_K(param_count+1)); 
        param_count = param_count+2;
    elseif model.B(j) == -1
            tmp_ODE_B(j) = tmp_ODE_B(j) - ...
            sym_K(param_count)*(sym_X(j))/((sym_X(j))+sym_K(param_count+1));             
            param_count = param_count+2;
    end
    % add the terms for basal regulations
    if model.C(j) == -1
        tmp_ODE_A(j) = tmp_ODE_A(j) - ...
            sym_K(param_count)*d*(sym_X(j))/((sym_X(j))+sym_K(param_count+1)); 
            param_count = param_count+2; 
    elseif model.C(j) == 1
        tmp_ODE_A(j) = tmp_ODE_A(j) + ...
            sym_K(param_count)*d*(1-sym_X(j))/((1-sym_X(j))+sym_K(param_count+1)); % need to make sure this is correct
            param_count = param_count+2;
    end
end
    
    % add the terms for the internode regulations
for j=1:n_nodes
    for k=1:n_nodes
        if model.A(j,k) == 1
            tmp_ODE_A(j) = tmp_ODE_A(j) + ...
            sym_K(param_count)*sym_X(k)*(1-sym_X(j))/((1-sym_X(j))+sym_K(param_count+1)); % need to make sure this is correct
            param_count = param_count+2; 

        elseif model.A(j,k) == -1
            tmp_ODE_A(j) = tmp_ODE_A(j) - ...
            sym_K(param_count)*sym_X(k)*(sym_X(j))/((sym_X(j))+sym_K(param_count+1)); % need to make sure this is correct
            param_count = param_count+2; 
        end            
    end   
    
end
sym_vec = [t; d; sym_K; sym_X];
model.param_count=n_parameters;
model.sym_ODE = tmp_ODE_A + s* tmp_ODE_B; 
% neet to convert the symbolic ODEs to the necessary functions
args = {t,sym_X,sym_K,d};
model.tmp_ODE_A = tmp_ODE_A;
tmp_ODE_AA = matlabFunction(tmp_ODE_A,'Vars',args);
tmp_ODE_BB = matlabFunction(tmp_ODE_B,'Vars',args);
%model.ODE = @(t,x,k,d) tmp_ODE_AA(t,x,k,d)+input_profile(t)*tmp_ODE_BB(t,x,k,d);
%model.ODE = @(t,x,k,d) eval(subs(tmp_ODE_A,sym_vec,[t d k x']'))+input_profile(t)*eval(subs(tmp_ODE_B,sym_vec,[t d k x']'));  

Jacobian_A = matlabFunction(jacobian(tmp_ODE_A,sym_X),'Vars',args);
Jacobian_B = matlabFunction(jacobian(tmp_ODE_B,sym_X),'Vars',args);
model.Jacobian = @(t,x,k,d) Jacobian_A(t,x,k,d)+input_profile(t)*Jacobian_B(t,x,k,d);

% Find the derivative of the ODE with respect to the parameters. 
K_matrix_A = matlabFunction(jacobian(tmp_ODE_A,sym_K(parameter)),'Vars',args);
K_matrix_B = matlabFunction(jacobian(tmp_ODE_B,sym_K(parameter)),'Vars',args);

% sys ODEs only sees the first n_nodes of x
sys_ODEs = @(t,x,k,d) tmp_ODE_AA(t,x,k,d)+input_profile(t)*tmp_ODE_BB(t,x,k,d);
model.sys_ODEs = sys_ODEs;
% sens ODEs uses all x. 
sens_ODEs = @(t,x,k,d) model.Jacobian(t,x(1:n_nodes),k,d)*x(n_nodes+1:end)+ ...
    K_matrix_A(t,x(1:n_nodes),k,d)+input_profile(t)*K_matrix_B(t,x(1:n_nodes),k,d);
% construct the full ODEs
model.ODE = @(t,x,k,d) [sys_ODEs(t,x(1:n_nodes),k,d) ; sens_ODEs(t,x,k,d)]; 
%model.ODE = @(t,x,k,d) eval(subs(tmp_ODE_A,sym_vec,[t d k x']'))+input_profile(t)*eval(subs(tmp_ODE_B,sym_vec,[t d k x']'));  

end  
            
    