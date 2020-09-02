function OBJ = get_simp_OBJ0(kinetics_models,kinetics_data,all_params,kk)
    basal_param = 0.1;    

    OBJ = 0; 
    for i=1:length(kinetics_data)
        
        m = kinetics_models{i}; 
        data = kinetics_data{i};
        
        ODE = @(t,x)m.model.ODE(t,x,all_params,basal_param);
        JAC = @(t,x)m.model.Jacobian(t,x,all_params,basal_param);
        IC = 0.05*ones(m.model.n_nodes,1); IC(end)=0; 
        options = odeset('Jacobian',JAC);
         
        [~,yout] = ode23s(ODE,data.tt*60,IC,options);
        model_hog = yout(:,m.model.n_nodes);
        OBJ = OBJ + nanmean(((data.hogp-model_hog)./data.STDVHog).^2);
    drawnow;subplot(2,2,kk); hold on; plot(data.hogp); plot(model_hog); 
    end
    OBJ = OBJ/length(kinetics_data); 
end