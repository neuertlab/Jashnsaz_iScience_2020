function [HOGS] = get_predictions(model,par_chain,data)           

[Hog1pp] = @(params) get_simp_sol(model,params,data);

ii  = 1; 
for i=1:size(par_chain,1)
    pars = par_chain(i,:)'; 
    try
        HOGS(ii,:) = Hog1pp(pars);
        ii = ii + 1; 
    catch
        disp('ode fail ...') % to avoid potential ODE solver failure
    end
end

end