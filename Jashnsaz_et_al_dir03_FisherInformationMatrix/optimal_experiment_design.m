function optimal_experiment_design()
parpool(4); 
spmd
    warning('off','all')
end

% models
load Models_TrueParams_simData/Models
model=Model{3}; % the true model

% true parameters 
load Models_TrueParams_simData/FP_OBJ_FIM.mat
true_params = FP_OBJ_FIM.best_pars;  clear FP_OBJ_FIM; 

% simulated data
load Models_TrueParams_simData/sim_data 
data_pool = [1:36]; 

%% calculate FIM matrices for individual pool of 36 datasets or load them
free_params = [1:model.n_params]; n_params = model.n_params; 
if exist('FIMs.mat')~=2
    log_e_scale=1; % to calculate FIM in log parameter space
    parfor i=1:length(data_pool)
        data=sim_data{i}; 
        [~,FIM_Matrix,~] = get_FIM(model,true_params,free_params,data,log_e_scale, 0);
        FIMs(:,:,i) = FIM_Matrix;         
    end
    save('FIMs','FIMs'); 
else
    load FIMs;
end
disp(' ')
disp(' ')

%% experiment designs
kinetics_types = {'STEPS', 'DIVERSE07M','DOPTIMAL'};
DATA_IDs = [[1:6];[6:6:36]]; % Figure 1H

cols = {'353848', '9baec8', 'A593E0'};
for i=1:length(cols)
    cmap(i,:)=hex2rgb(cols{i}); 
end
select_expmnts=[1 2 3]; 

figure(1); clf; set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 8 4], 'PaperUnits', 'centimeters', 'PaperSize', [8 4]); 
dx=.2; dy=.2; 
subplotHJ(1,1,1,dy,dx); hold on; axis on; box on; 
FaceCol = [1 1 1]*.75; scFaceAlpha=.15; sz=35; 

%% optimal experiment design; 36 choose n

I0=eye(n_params); 

for n=1:6 % for number of datasets 1 to 6
    EXPMNTS_ID  = []; 
    proposed_experiments = nchoosek(data_pool,n); % propose all combinatorials of n choose 36 datasets
    tic
    parfor i=1:size(proposed_experiments,1)
        FIM = nansum(FIMs(:,:,[proposed_experiments(i,:)]),3); % sum FIMs for proposed experiments
        invFIM = I0/FIM; % inverse sum of FIM 
        DETinvFIMs(i) = det(invFIM); % get determinant of FIM^{-1}
    end
    disp(['calculation time: ',num2str(toc), ' seconds.'])
    DETinvFIMs=abs(DETinvFIMs)/(log(10))^(2*n_params); % convert from  natural log to log10
    [~,INDX]=sort(DETinvFIMs,'ascend'); 
    RESULTS.proposed_experiments{n} = proposed_experiments; % datasets IDs in all combinatorial proposed experiments
    RESULTS.DETinvFIMs{n} = DETinvFIMs; % determinants of FIM^{-1} for all experiments
    RESULTS.INDX{n} = INDX; % sorted indexes (lowest to highest determinant value) for all experiments
    
%     DETinvFIMs = RESULTS.DETinvFIMs{n}; 
%     INDX = RESULTS.INDX{n};

    DOPTIMAL_EXPMNT_INDX = INDX(1); % FIM D_Optimal Experiment (Experiment that minimizes determinant of FIM^{-1})

    for i=1:size(DATA_IDs)
        EXPMNTS_ID(i) = find(ismember(proposed_experiments,DATA_IDs(i,1:n),'rows')==1); 
    end
    EXPMNTS_ID = [EXPMNTS_ID DOPTIMAL_EXPMNT_INDX]; 
    disp(['D_Optimal Experiment IDs of n=', num2str(n), ', among all 36_choose_n (=', num2str(length(INDX)), ') possible experiments:']) 
    disp(proposed_experiments(DOPTIMAL_EXPMNT_INDX,:))
    
    % plot determinant FIM^{-1} versus number of data sets 
%     scatter(n + .1*randn(length(DETinvFIMs),1), DETinvFIMs, 2, 'o', 'MarkerEdgeColor','none', 'MarkerFaceColor', FaceCol, 'MarkerFaceAlpha', scFaceAlpha); 
    cmapp=cmap([select_expmnts],:); 
    ii=1;
    for i=EXPMNTS_ID(select_expmnts)
        scatter(n, DETinvFIMs(i), 22, 'o','MarkerEdgeColor','none','MarkerFaceColor', cmapp(ii,:));
        ii=ii+1; 
    end
    xlim([.5 6.5]); xticks([1:6]);
    set(gca,'YScale', 'log'); yticks([1e0 1e50 1e100]);
    xlabel('Number of train data, n'); ylabel('Det{FIM^{-1}}')
    legend(kinetics_types)

    set(findall(gcf,'-property','FontSize'),'FontSize',7, 'defaultTextFontSize',7, 'FontName', 'Helvetica');

    figname = ['DETinvFIM_vs_n']; 
    print(figname,'-depsc', '-r600'); 
    
end
save('RESULTS','RESULTS'); 

poolobj = gcp('nocreate');
delete(poolobj);

end


