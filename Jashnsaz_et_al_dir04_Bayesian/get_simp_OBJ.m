function OBJ = get_simp_OBJ(kinetics_models,kinetics_data,params,free_params,all_params)
    all_params(free_params) = params; % pass on free model params, and rest the rest (if any) to the true parameters
    basal_param = 0.1; % basal paramters

    OBJ = 0; 
    % loop through all the training data set
    for i=1:length(kinetics_data)
        m = kinetics_models{i}; % model
        data = kinetics_data{i}; % data       
        ODE = @(t,x)m.model.ODE(t,x,all_params,basal_param); % model ODE
        JAC = @(t,x)m.model.Jacobian(t,x,all_params,basal_param); % model Jacobian
        IC = 0.05*ones(m.model.n_nodes,1); IC(end)=0; % initial condition
        options = odeset('Jacobian',JAC); %'AbsTol', 1e-6, 'RelTol', 1e-3,  'NormControl','on'  
        clear t_obj; t_obj=tic; 
        [~,yout] = ode23s(ODE,data.tt*60,IC,options);        
        if toc(t_obj)<0.1 % to avoid potential ode23s solver stock
            model_hog = yout(:,m.model.n_nodes);    
            try
                OBJ = OBJ + nanmean(((data.hogp-model_hog)./data.STDVHog).^2);
            catch
                OBJ = OBJ + inf; 
            end
        else
            OBJ = OBJ + inf; 
        end
    end
    OBJ = OBJ/length(kinetics_data); 
end