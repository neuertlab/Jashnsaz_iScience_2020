function [simData] = get_simData()

dir_name='simData1'; 
mkdir(dir_name);

% load models and their best parameters to generate synthetic data
load Models_and_TrueParams/Models
model_index = 3; % the true model
model=Model{model_index}; 

load Models_and_TrueParams/FP_OBJ_FIM.mat % get true parameters
params = FP_OBJ_FIM.best_pars; 

% set total treatment duration (TTT) and treatment ON duration (TT)
TTT = 50; TT=25; 

% pick functions from roots and polynomials
roots = [10 5 3 2 1.5]; 
polynomials = [1:14]; 

% Final NaCl molarities (M)
NaCls = [1e-3 2e-3 5e-3 1e-2 2e-2 5e-2 [.1 .2 .3 .4 .5 .6 .7 .8 .9 1] 2 5 10 20]; 

n_polynomials = length(polynomials); 
n_roots = length(roots);
n_NaCls = length(NaCls); 
all_funs = 1+n_polynomials+n_roots; 
tot_inputs = all_funs*n_NaCls+1; 

fxngrps = 1; 
% step function
sympref('HeavisideAtOrigin',1);
step_fun = @(t, maxNaCl) maxNaCl * heaviside(t - .001);

fxngrps = fxngrps + 1; 
% root functions
for i=1:n_roots 
    root_fun{i} = @(t, maxNaCl) min([(maxNaCl.*t.^(1/roots(i))/(60*TT)^(1/roots(i))),maxNaCl]); 
    fxngrps = fxngrps + 1; 
end
lineramp = fxngrps; 

% polynomial functions
for i=polynomials
    poly_fun{i} = @(t, maxNaCl) min([(maxNaCl.*t.^i/(60*TT)^i),maxNaCl]);
    fxngrps = fxngrps + 1; 
end

%control NaCl=constant; 
contrl_fun = @(t, maxNaCl) maxNaCl + (t - t);

% build the input profiles for all kinetic types and final concentrations
fun_n = 1;
for steps=NaCls
    salt_funcs{fun_n} = @(t) step_fun(t, steps); 
    fun_n = fun_n +1; 
end

for rootfuns = 1:n_roots
    for ramps=NaCls
        salt_funcs{fun_n} = @(t) root_fun{rootfuns}(t, ramps); 
        fun_n = fun_n +1; 
    end
end

for polys = polynomials
    for ramps=NaCls
        salt_funcs{fun_n} = @(t) poly_fun{polys}(t, ramps); 
        fun_n = fun_n +1; 
    end
end

%control NaCl=0; 
salt_funcs{fun_n} = @(t) contrl_fun(t, 0); 


% names
saltfunctions{1} = '$$ steps $$'; ii=2; 
for i=roots 
    saltfunctions{ii} = ['$$ t^{1/', num2str(i), '} $$']; ii = ii +1; 
end

for i=polynomials 
    saltfunctions{ii} = ['$$ t^{', num2str(i), '} $$']; ii = ii +1; 
end


figure(2); set(gcf,'defaultLineLineWidth',1); 
set(gcf, 'Units', 'points', 'Position', [0 0 1600 1000], 'PaperUnits', 'points', 'PaperSize', [1600, 1000])
lw = 1;   
title_x_loc = 42; 

%% get simulated data for all input profiles
n_synthetic_replicates = 1; 
for sdata=1:n_synthetic_replicates
    simData(sdata).saltfunctions = saltfunctions; 
    simData(sdata).sim_data_params = params;  

    for i = 1:tot_inputs
    
    data{i}.tt=[0:TTT];
    data{i}.Salt = salt_funcs{i};
    model = Get_ODE(model,data{i}.Salt);                           % add ODE to model object for a given Salt input
    model.IC = 0.05*ones(model.n_nodes,1); model.IC(end)=0;        % set initial condition 

    % Define parameters.     
    kinetic_params = params; 
    basal_param=0.1;
    % Specify ODE
    ODE = @(t,x)model.ODE(t,x,kinetic_params,basal_param); 
    JAC = @(t,x)model.Jacobian(t,x,kinetic_params,basal_param);
    options = odeset('Jacobian',JAC);
     % Integrate ODE. 
    [~,yout] = ode23s(ODE,data{i}.tt*60,model.IC);

    data{i}.MeanHog = yout(:,model.n_nodes);                       % get the hog output - the "bottom" node. 
    data{i}.cumsumMeanHog = cumsum(data{i}.MeanHog);  

    %generate single cell trajs by adding Gaussian noise to the mean.
    noise_level = 0.02; 
    BiolRep = 5; 
    single_cell_trajs = 10; 

    for br=1:BiolRep
        mean_var=1+0*randn(1,1);
        for sc=1:single_cell_trajs
            data{i}.scSTDVHog(:,sc) = noise_level*randn(size(yout,1),1); %mean(yout(:,model.n_nodes))
            data{i}.scHogp(:,sc) = mean_var*yout(:,model.n_nodes) + data{i}.scSTDVHog(:,sc);  
            data{i}.cumsum_scHogp(:,sc) = cumsum(data{i}.scHogp(:,sc));
        end
        data{i}.br_hogp(:,br) = mean(data{i}.scHogp')';
        data{i}.br_STDVHog(:,br) = std(data{i}.scHogp')';
        data{i}.cumsum_br_hogp(:,br) = mean(data{i}.cumsum_scHogp')';
        data{i}.cumsum_br_STDVHog(:,br) = std(data{i}.cumsum_scHogp')';
    end

    data{i}.hogp = mean(data{i}.br_hogp')';
    data{i}.STDVHog = 2*std(data{i}.br_STDVHog')';

    data{i}.cumsumhogp = mean(data{i}.cumsum_br_hogp')';
    data{i}.cumsumSTDVHog = std(data{i}.cumsum_br_STDVHog')';

    % calculate salt function for given time points.
    j=1;
    salt_out=zeros(1,length(data{i}.tt));
    for t = data{i}.tt
        salt_out(j) = data{i}.Salt(t*60);
        j=j+1;
    end 
    data{i}.salt_out=salt_out;

    end
    simData(sdata).data = data; 

% 6x6 sim_data 
    pick_NaCls = [.05 .1 .2 .3 .5 .7]; 
    pick_functions = [1 5 7 8 11 13]; 

    iii = 1; 
    ii = 1; 
    for i = 1:all_funs
    for j=1:n_NaCls
        if (ismember(NaCls(j), pick_NaCls)==1 & ismember(i,pick_functions)==1)
            sim_data{ii} = data{iii}; 
            ii = ii + 1;
        end
        iii = iii + 1; 
    end
    end

    for i=1:length(pick_functions) 
        sim_data_titles{i} = saltfunctions{pick_functions(i)};  
    end
    % control (0 input over time)
    sim_data{ii} = data{tot_inputs};
    
    simData(sdata).sim_data=sim_data;  
    simData(sdata).simdata_maxNaCls = pick_NaCls; 
    simData(sdata).simdata_funcs = pick_functions; 

% sim_test_data1
    pick_NaCls = [.9 1];
    pick_functions = [1 5 7 8 11 13]; 

    iii = 1; 
    ii = 1; 
    for i = 1:all_funs
        for j=1:n_NaCls
            if (ismember(NaCls(j), pick_NaCls)==1 & ismember(i,pick_functions)==1)
                sim_test_data1{ii} = data{iii}; 
                ii = ii + 1;
            end
            iii = iii + 1; 
        end
    end
    simData(sdata).sim_test_data1=sim_test_data1; 
    simData(sdata).testdata1_maxNaCls = pick_NaCls; 
    simData(sdata).testdata1_funcs = pick_functions; 

% sim_test_data2
    pick_NaCls = [0.05 .1 .2 .3 .5 .7];
    pick_functions = [15];

    iii = 1; 
    ii = 1; 
    for i = 1:all_funs
        for j=1:n_NaCls
            if (ismember(NaCls(j), pick_NaCls)==1 & ismember(i,pick_functions)==1)
                sim_test_data2{ii} = data{iii}; 
                ii = ii + 1;
            end
            iii = iii + 1; 
        end
    end
    simData(sdata).sim_test_data2=sim_test_data2;
    simData(sdata).testdata2_maxNaCls = pick_NaCls; 
    simData(sdata).testdata2_funcs = pick_functions; 
 
    for i=1:length(pick_functions) 
        sim_test_data2_titles{i} = saltfunctions{pick_functions(i)};  
    end

%     clear data; 
% end

save([dir_name,'/simData'],'simData'); 
disp(['simData(' num2str(sdata), ') is simulated from Model', num2str(model_index),'.']);
    
%% ploting

    % plot salts (inputs)
    if sdata == 1
        figure(1); dx=0.05; dy=0.05;
        set(gcf,'defaultLineLineWidth',.5); 
        set(gcf, 'Units', 'centimeters', 'Position', [0 0 21 21], 'PaperUnits', 'centimeters', 'PaperSize', [21 21])
        set(gca, 'FontName', 'Helvetica'); 
        set(gca,'defaultAxesTickLabelInterpreter','latex');
        cmap = winter(6); %winter bone gray copper
        sim_data_titles = {'steps $ (t^{0}) $', 'root2 $ (t^{1/2}) $', 'linears $ (t^{1}) $', 'quadratics $ (t^{2}) $', 'quintics $ (t^{5}) $', 'heptics $ (t^{7}) $'}; 

        pickdata = 1; 
        for i = 1:6
            subplotHJ(3,3,i,dy,dx); grid on
            for j=1:6   
                hold on
                plot(sim_data{pickdata}.tt,sim_data{pickdata}.salt_out, 'LineWidth', 3, 'color', cmap(j,:));
                pickdata = pickdata + 1; 
                box off   
            end
            txt = text(title_x_loc,.95, [sim_data_titles{i}], 'Interpreter','latex'); txt.HorizontalAlignment = 'center'; 
            xlim([-.5,inf]); ylim([0,1]); xticks([0:10:50])          
            ylabel('NaCl (M)'); xlabel('time (min)'); box on
        end    

        % sim_test_data1 
%         dy = 1.5*dy; 
        sim_test_data1_titles = {'$ 0.9M $', '$ 1.0M $'}; 
        pickdata = 1; 
        for i = 1:2
            subplotHJ(3,3,6+i,dy,dx); %grid on
            for j=1:6   
                hold on
                if i==1 
                    pickdata = 2*j-1; 
                    txt = text(35,.2+j*0.05, sim_data_titles{j}, 'Interpreter','latex'); 
                    txt.HorizontalAlignment = 'left'; txt.Color = cmap(j,:);
                else
                    pickdata = 2*j; 
                end
                plot(sim_test_data1{pickdata}.tt,sim_test_data1{pickdata}.salt_out, 'LineWidth', 3, 'color', cmap(j,:));
                box off   
            end
            txt = text(title_x_loc,.95, [sim_test_data1_titles{i}], 'Interpreter','latex', 'Fontsize', 5); txt.HorizontalAlignment = 'center'; 
            xlim([-.5,inf]); ylim([0,1]); xticks([0:10:50])          
            ylabel('NaCl (M)'); xlabel('time (min)'); box on
        end    

        % sim_test_data2 
        sim_test_data2_titles = {'nonics $ (t^{9}) $'}; 
        pickdata = 1; 
        subplotHJ(3,3,9,dy,dx); grid on

        for i = 1:1
            for j=1:6   
                hold on
                plot(sim_test_data2{pickdata}.tt,sim_test_data2{pickdata}.salt_out, 'LineWidth', 3, 'color', cmap(j,:)-(i-1)*.05*cmap(j,:));
                pickdata = pickdata + 1; 
                box off   
            end    
        end 

        txt = text(title_x_loc,.95, [sim_test_data2_titles], 'Interpreter','latex'); 
        txt.HorizontalAlignment = 'center'; txt.VerticalAlignment = 'top'; 
        xlim([-.5,inf]); ylim([0,1]); xticks([0:10:50])          
        ylabel('NaCl (M)'); xlabel('time (min)'); box on

        set(findall(gcf,'-property','FontSize'),'FontSize',7.5, 'defaultTextFontSize',7.5, 'FontName', 'Helvetica')
        
        print([dir_name, '/inputs'],'-dpng');
        saveas(gcf,[dir_name, '/inputs'], 'epsc') 
    end
%% plot Hog1
 hh=figure(2); dx=0.03; dy=0.05;

 p1 = winter(100); 
 p2 = summer(100); 
 cmap2 = [parula(36); [1 0 0]; p1(100-11:100,:); p2(1:6,:)];

 pns = [1:3 5+(1:3) 10+(1:3)];  
 sim_NaCls = [0.05 .1 .2 .3 .5 .7];
    ccc=1; 
    pickdata = 1; 
    for i = 1:6
        subplotHJ(3,5,pns(i),dy,dx);
        for j=1:6   
            hold on
            if i==1 & sdata==1
                txt = text(40,.2+j*0.05, [num2str(sim_NaCls(j), '%4.2f'), 'M'], 'Interpreter','latex'); 
                txt.HorizontalAlignment = 'left'; txt.Color = cmap(j,:);
            end
            shadedErrorBar(sim_data{pickdata}.tt, sim_data{pickdata}.hogp, sim_data{pickdata}.STDVHog, {'LineWidth',lw, 'color', cmap(j,:)}, 0.15);
            pickdata = pickdata + 1; ccc = ccc + 1; 
        end
        if i==3; title(['BR ', num2str(sdata)]); end
%         if i==1;             
%             shadedErrorBar(sim_data{31}.tt, sim_data{31}.hogp, sim_data{31}.STDVHog, {'LineWidth',lw, 'color', cmap(end,:)}, 0.15);
%         end

        if sdata == 1; 
            txt = text(title_x_loc,.6, [sim_data_titles{i}], 'Interpreter','latex'); txt.HorizontalAlignment = 'center'; 
        end
        xlim([-.5,inf]); ylim([0,.65]); xticks([0:10:50])          
        ylabel('Hog1pp'); xlabel('time (min)'); box on
    end    
    
    ccc = ccc + 1; 
    
    sim_test_data1_titles = {'0.9M', '1.0M'}; 
    pickdata = 1; 
    for i = 1:2
        subplotHJ(3,5,pns(6+i),dy,dx)
        for j=1:6   
            hold on
            if i==1; 
                pickdata = 2*j-1; 
                if sdata == 1; 
                    txt = text(35,.2+j*0.05, sim_data_titles{j}, 'Interpreter','latex'); 
                    txt.HorizontalAlignment = 'left'; txt.Color = cmap(j,:);
                end
            else
                pickdata = 2*j; 
            end
            shadedErrorBar(sim_test_data1{pickdata}.tt, sim_test_data1{pickdata}.hogp, sim_test_data1{pickdata}.STDVHog, {'LineWidth',lw, 'color', cmap(j,:)}, 0.15);
            ccc = ccc + 1; 
        end
        if sdata == 1; 
            txt = text(title_x_loc,.6, [sim_test_data1_titles{i}], 'Interpreter','latex'); txt.HorizontalAlignment = 'center'; 
        end
        xlim([-.5,inf]); ylim([0,.65]); xticks([0:10:50])
        ylabel('Hog1pp'); xlabel('time (min)', 'Fontsize', 6); box on
    end    
    
    sim_test_data2_titles = {'nonics $ (t^{9}) $'}; 
    pickdata = 1; 
    subplotHJ(3,5,pns(end),dy,dx)
    tst_NaCls = [0.05 .1 .2 .3 .5 .7];
    for i = 1:1
        for j=1:6   
            hold on      
            if i==1 
                if sdata == 1; 
                    txt = text(40,.2+j*0.05, [num2str(tst_NaCls(j),'%4.2f'), 'M'], 'Interpreter','latex'); 
                    txt.HorizontalAlignment = 'left'; txt.Color = cmap(j,:);
                end
            end
            shadedErrorBar(sim_test_data2{pickdata}.tt, sim_test_data2{pickdata}.hogp, sim_test_data2{pickdata}.STDVHog, {'LineWidth',lw, 'color', cmap(j,:)}, 0.15);           
            pickdata = pickdata + 1; ccc = ccc + 1; 
        end   
    end 
    if sdata == 1; 
        txt = text(title_x_loc,.6, [sim_test_data2_titles], 'Interpreter','latex'); 
        txt.HorizontalAlignment = 'center'; txt.VerticalAlignment = 'top';
    end
    xlim([-.5,inf]); ylim([0,.65]); xticks([0:10:50])          
    ylabel('Hog1pp'); xlabel('time (min)'); box on
    
    %% plot noise
    pns = []; pns = [4 9 14]; 
    
    count = 1; 
        for i=1:length(sim_data)
            all_sim_data{count} = sim_data{count}; 
            count = count + 1; 
        end
        for i=1:length(sim_test_data1)
            all_sim_data{count} = sim_test_data1{i}; 
            count = count + 1; 
        end
        for i=1:length(sim_test_data2)
            all_sim_data{count} = sim_test_data2{i}; 
            count = count + 1; 
        end

        
        subplotHJ(3,5,pns(1),dy,dx); 
        for i=1:length(all_sim_data)
            hold on
            meannoise(i) = mean(all_sim_data{i}.STDVHog); 
            plot(i, meannoise(i),'--*', 'Color', cmap2(i,:))
        end
        plot(meannoise,':k'); 
        
        xlim([0 inf]); ylim([7e-3 11e-3])
        xlabel('sim data'); ylabel('noise')
        box on; xticks([0:10:90]); %title('D1 - D53'); 

        subplotHJ(3,5,pns(2),dy,dx); 
        for i=1:length(all_sim_data)
            hold on
            meannoise(i) = mean(all_sim_data{i}.hogp)./mean(all_sim_data{i}.STDVHog); 
            plot(i, meannoise(i),'--*', 'Color', cmap2(i,:))
        end
        if sdata==1; plot(meannoise,':k'); end
        ylim([0 30]); xlim([0 inf]); 
        xlabel('sim data'); ylabel('signal to noise')

        box on; xticks([0:10:90]); 
        
        sp = subplotHJ(3,5,pns(3),dy,dx); cla(sp); 
        for i=1:length(all_sim_data)
            hold on
            plot(all_sim_data{i}.hogp, all_sim_data{i}.STDVHog, 'o', 'Color', cmap2(i,:))
        end
        set(gca, 'XScale', 'log' ); 
        xlim([1e-6 1]); ylim([0 .025]);
        
        xlabel('Hog1pp'); ylabel('noise'); box on
         
    %% get dose response 
    pns = []; pns = [5 10 15]; 
    dose_response = NaN(all_funs,n_NaCls); 
    max_activation_time = NaN(all_funs,n_NaCls); 

    for i = 1:all_funs
        for j=1:n_NaCls
            pickdata = (i-1)*n_NaCls+j; 
            dose_response(i,j) = max(data{pickdata}.hogp); 
            cum_active(i,j) = mean(data{pickdata}.hogp); 

            max_activation_time(i,j) = median(data{pickdata}.tt(find(data{pickdata}.hogp>.9*max(data{pickdata}.hogp)))); 
            tmaxHog = floor(max_activation_time(i,j)); 
            if tmaxHog == 0
                max_activation_time(i,j)= NaN; 
            elseif mean(data{pickdata}.hogp([tmaxHog]))<1.1*mean(data{pickdata}.hogp); 
                max_activation_time(i,j) = NaN; 
            end
        end   
    end  
%     cum_active = cum_active./max(max(cum_active));

    simData(sdata).dose_response = dose_response; 
    simData(sdata).cum_active = cum_active; 
    simData(sdata).tmax_activ = max_activation_time; 
    save([dir_name,'/simData'],'simData'); 

    % plot dose and activation time
    cmap3=hot(20); 
    lineramp = 7; 
        
    subplotHJ(3,5,pns(1),dy,dx); 
    for i=1:all_funs
        hold on; if i==lineramp; col = [0 .9 .7]; else col=cmap3(i,:); end
        plot(NaCls, dose_response(i, :), ':*', 'LineWidth', 1, 'Color', col)
    end

    ylim([0 1]); xlim([NaCls(1) NaCls(end)])
    xticklabels([]); 
    set(gca, 'XScale', 'log' ); box on
    ylabel('Hog1pp_{max}')   
    lgd=legend(saltfunctions,'Interpreter','latex','Color', 'w', 'Location','northwest'); legend('boxoff');
    
    subplotHJ(3,5,pns(2),dy,dx); 
    for i=1:all_funs
        hold on; if i==lineramp; col = [0 .9 .7]; else col=cmap3(i,:); end
        plot(NaCls(1:end-3), cum_active(i, 1:end-3), ':*', 'LineWidth', 1, 'Color', col)
    end

    ylim([0 inf]); xlim([NaCls(1) NaCls(end)])
    xticklabels([]); 
    set(gca, 'XScale', 'log' ); box on

    ylabel(' \int{Hog1pp}')

    
    subplotHJ(3,5,pns(3),dy,dx); 
    for i=1:all_funs
        hold on; if i==lineramp; col = [0 .9 .7]; else col=cmap3(i,:); end
        plot(NaCls(4:end-2), max_activation_time(i, 4:end-2), ':*', 'LineWidth', 1, 'Color', col)
    end

    ylim([0 50]); xlim([NaCls(1) NaCls(end)])
    set(gca, 'XScale', 'log' ); box on

    xlabel('final NaCl [M]'); ylabel('Hog1pp_{max} time (min)')

    set(findall(gcf,'-property','FontSize'),'FontSize',11, 'FontName', 'Helvetica'); 
    lgd.FontSize = 8;
    
    print([dir_name, '/Responses', num2str(sdata)],'-dpng');
    saveas(gcf,[dir_name, '/Responses', num2str(sdata)], 'epsc') 

    myframes(sdata) = getframe(hh);
    clear data; 
end
    implay(myframes, 1);
        
    VideoName = [dir_name, '/Responses'];
    v = VideoWriter(VideoName); 
    v.FrameRate = 0.5;
    open(v)   
    writeVideo(v,myframes);
    close(v);
end