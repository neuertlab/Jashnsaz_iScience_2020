function Bayesian(cluster_vs_local,free_parameters,kinetics_name,train_data_IDs,chain)
warning('off','all'); 

switch cluster_vs_local
    case 'local'
        myCluster = parcluster('local');
        n_cores = myCluster.NumWorkers; % call # of cores on the local machine 
    case 'cluster'
        n_cores = str2num(getenv('SLURM_JOB_CPUS_PER_NODE')); % call # of cores on the cluster (SLURM)   
end

parpool(n_cores); % start parallel pool
spmd
    warning('off','all')
end

% variables of Bayesian algorithm
LUB = 3; % parameters upper and lower bound in log10 space; uniform prior [-LUB, +LUB] initialization
nSamples = 5000; % number of parameter sets to be sampled at each iteration
flattenOBJ = .01*nSamples; % flatten 1% of parameter sets OBJs to max of the 1%
iMax = 1000; % number of iterations

% models
load Models_TrueParams_simData/Models
model=Model{3}; % the true model

% true params 
load Models_TrueParams_simData/FP_OBJ_FIM; 
true_params = FP_OBJ_FIM.best_pars; clear FP_OBJ_FIM; 
n_pars  =  length(free_parameters); % # of free model parameters

% simulated data
load Models_TrueParams_simData/sim_data

% randomize the seed for rng (important to initiate different chains)
rng('shuffle'); 
seeds = randi(2^32,1e4,1); 
rng(seeds(randi(1e4,1))); 

%% get kinetics, define roots, make directories, define OBJ function 
disp([kinetics_name, ' kinetics, chain',num2str(chain), ' started  ...']); 
dirName =  [kinetics_name, '/','chain',num2str(chain),'/']; mkdir(dirName); disp(' ')

% training data and their models ODEs
for i=1:length(train_data_IDs)
    traindata{i} = sim_data{train_data_IDs(i)};
    traindata_models{i}.model = Get_ODE(model,traindata{i}.Salt);
end
BayesResults.train_data_IDs=train_data_IDs; 

% define the objective function for model and training data
kinetics_models = traindata_models; kinetics_data = traindata; 
[OBJ] = @(param) get_simp_OBJ(kinetics_models,kinetics_data,10.^param,free_parameters,true_params);

%% start the Bayesian 
    
    InitialPopulation = -LUB + 2*LUB*rand(nSamples,n_pars); % initialize parameters with uniform prior
    
    mu_old = mean(InitialPopulation); % mean parameters
    cov_old = cov(InitialPopulation); % covariance parameters
    eta = 0.5;  % fit update inertia
    eta2 = .01; % coefficient for the identity matrix to be added to the covariance matrix
    
    full_posterior = []; mu_s = []; cov_s =[]; 
    for i = 1:iMax
        ti=tic; 
        disp(['i = ', num2str(i), ' started.'])
        if i == 1; par_chain = InitialPopulation; end
        full_posterior = [full_posterior; par_chain];

        % solve the model and calculate OBJ for each parameter set
        OBJs = inf(1,nSamples); 
        parfor j=1:nSamples
%             if mod(j,50)==0; disp(['  ... params set ', num2str(j)]); end
            OBJs(j) = OBJ(par_chain(j,:)); 
        end     
        discards(i) = length(find(isinf(OBJs)==1)); % quantify the # parameter set that ODE solver was unsuccesful
                
        [sOBJs, INXs] = sort(OBJs, 'ascend'); % sort OBJs        
        % flatten the lowest 1% of OBJs to the highest OBJ among them
        OBJs(INXs(1:flattenOBJ)) = OBJs(INXs(flattenOBJ));
        nL=exp(-OBJs)/nansum(exp(-OBJs));  % likelihood of parameter sets

        % weighted mean and covariance based on the likelihood values
        mu = (nL*par_chain);  % weighted mean of parameters
        ss = nL.*(par_chain-mu)'*(par_chain-mu);  % weighted coveriance of parameters
        
        % add inertia terms to mu and cov
        mu = eta*mu+(1-eta)*mu_old;
		ss = eta*ss+(1-eta)*cov_old;
		if i<iMax/2
			ss = ss + eta2*eye(n_pars);
		end
		
		mu_old = mu;
		cov_old = ss;
        
        w_mean_pars = mu; w_cov_pars = .5*(ss+ss');         
%         [fig] = plot_MVN(par_chain,mu,w_cov_pars,free_parameters,i,dirName);

        % re-sample parameter sets      
        par_chain = mvnrnd(w_mean_pars', w_cov_pars, nSamples); % sample parameter sets using mvnrnd
%         [fig] = plot_MVN(par_chain,mu,ss,free_parameters,i,dirName);

        par_chain(par_chain<-LUB)=-LUB; par_chain(par_chain>LUB)=LUB; % bound the sampled parameter sets to [-LUB, +LUB] 
        mu_s = [mu_s; w_mean_pars]; cov_s =[cov_s; w_cov_pars]; % collect the mu and cov at each iteration to save
        chains_OBJs(i)=sOBJs(1); % collect the best OBJ at each iteration
        disp(['best OBJ value is ', num2str(sOBJs(1))]);
        disp(['time = ', num2str(toc(ti)), 's | i = ', num2str(i), ' finished.'])
        disp(' ')  
        
    end
    BayesResults.InitialPopulation=InitialPopulation; 
    BayesResults.full_posterior=full_posterior; 
    BayesResults.w_mean_pars=w_mean_pars; 
    BayesResults.w_cov_pars=w_cov_pars; 
    BayesResults.chains_OBJs = chains_OBJs; 
    BayesResults.discards = discards; 
    BayesResults.mu_s = mu_s; 
    BayesResults.cov_s =cov_s; 
    save([dirName,'/BayesResults'], 'BayesResults');  
    disp(' ')      

    %% draw samples from the posterior to get  predictions
    disp(['draw predictions for 100 samples from the the posterior ... ']); tic
    mean_post = w_mean_pars; 
    cov_post = w_cov_pars; 
    par_samples = mvnrnd(mean_post, cov_post, 100); % sample parameter sets using mvnrnd

    all_params = repmat(true_params',100,1); 
    all_params(:,free_parameters) = 10.^par_samples;

    % predictions 
    parfor testdata=1:length(sim_data)       
        try
        [HOGS] = get_predictions(model,all_params,sim_data{testdata});      
        catch
            disp([' ... ODE with testdata',num2str(testdata), ' u(t) not solvable ...']); 
        end
        predictions{testdata}.HOGS = HOGS; 
        disp(['   response predictions on simdata ', num2str(testdata), ' is done.'])      
    end    
    disp('');
    disp('');

    BayesResults.predictions = predictions; clear predictions; 

    save([dirName,'/BayesResults'], 'BayesResults');  
    disp([' time = ', num2str(toc/60), 'min, predictions done.'])
    disp(' ')   

poolobj = gcp('nocreate');
delete(poolobj);

end


