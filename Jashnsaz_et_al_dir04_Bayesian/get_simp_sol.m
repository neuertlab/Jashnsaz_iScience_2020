function [Hog1pp] = get_simp_sol(model,params,data)

    model.IC = 0.05*ones(model.n_nodes,1); model.IC(end)=0;                
    basal_param = 0.1;   
    
    model = Get_ODE(model,data.Salt);
    ODE = @(t,x)model.ODE(t,x,params,basal_param);
    JAC = @(t,x)model.Jacobian(t,x,params,basal_param);
    options = odeset('Jacobian',JAC);
    
    try
    [~,yout] = ode23s(ODE,data.tt*60,model.IC,options);
    Hog1pp = yout(:,model.n_nodes);
    catch
        disp('ode23s failed ... ')
    end
end
