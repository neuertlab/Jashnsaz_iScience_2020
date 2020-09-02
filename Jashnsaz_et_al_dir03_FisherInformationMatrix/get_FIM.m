function [FIM,sFIM,S]= get_FIM(model,theta,free_parameters,data,log_e_scale, FIM0)
    % Get the Fisher Information Matrix for a particular model,
    % about parameters specified by theta. 
    % free parameters should be a vector of indices of theta that are 
    % "free" in the model. The variance should be of size Nt.
    % for now, assuming only observable is the last node. 
    % log_e_scale refers to the log-FIM (i.e. df/dlogtheta)^2
   
    % FOR each free parameter, compute the FIM:
    S = zeros(length(data.tt),length(free_parameters));
    for i=1:length(free_parameters)
        model = Get_Sens_ODE(model,data.Salt,free_parameters(i));
        x0 = zeros(model.n_nodes*2,1); x0(1:model.n_nodes-1,1) = 0.05;         
        basal_param = 0.1;
        ODE = @(t,x)model.ODE(t,x,theta,basal_param);
        [~,yout] = ode23s(ODE,data.tt*60,x0);
        S(:,i) = yout(:,end); % sensitivity matrix
    end
    % make the variance matrix
    COV = inv(diag(data.STDVHog.^2));
    % compute the FIM
    FIM = zeros(length(free_parameters),length(free_parameters));
    for i=1:length(free_parameters)
        for j=1:length(free_parameters)
            if log_e_scale
                FIM(i,j) = theta(free_parameters(i))*theta(free_parameters(j))*...
                    S(:,i)'*COV*S(:,j);
            else
                FIM(i,j) = S(:,i)'*COV*S(:,j);
            end            
        end 
    end
    sFIM = FIM + FIM0; 
    
    for i=1:size(sFIM,1)
        for j=1:size(sFIM,1)        
            if i>j
                sFIM(i,j) = sFIM(j,i);   
            end        
        end
    end 
    
end
