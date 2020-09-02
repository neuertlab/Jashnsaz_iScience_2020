function visualize_models()
clc
clear all
mkdir('ModelsGraphs'); 

load Models
models = [1:length(Model)]; 

h(1) = figure(1); set(gcf, 'Units', 'centimeter', 'Position', [0, 0, 40, 7], 'PaperUnits', 'centimeter', 'PaperSize', [40, 7]); 
set(gcf,'defaultLineLineWidth',1); set(0,'defaultAxesFontSize',7); set(gca, 'FontName', 'hevetica');

h(2) = figure(2); set(gcf, 'Units', 'centimeter', 'Position', [0, 0, 40, 40], 'PaperUnits', 'centimeter', 'PaperSize', [40, 40]); 
set(gcf,'defaultLineLineWidth',1); set(0,'defaultAxesFontSize',7); set(gca, 'FontName', 'hevetica');

count = 1;
for m=models
    figure(1); dx=0.015; dy=0.1; subplotHJ(1,5,count,dy,dx);
    network_plot(Model{m},m)
    m_fig=['Model', num2str(m), '.png'];
    image(imread(['ModelsGraphs/' m_fig]));
    title({['Model ', num2str(m), ' | N_{Params}=', num2str(Model{m}.n_params)]});
    axis normal
    axis off

    figure(2);  dx=0.06; dy=0.08; subplotHJ(1,1,1,dy,dx);
    network_plot(Model{m},m)
    m_fig=['Model', num2str(m), '.png'];
    image(imread(['ModelsGraphs/' m_fig]));
    title({['Model ', num2str(m), ' | N_{Params}=', num2str(Model{m}.n_params)]});
    axis normal
    axis off
    print(['ModelsGraphs/', m_fig(1:end-4)],'-dpng');

    count  = count + 1; 
end

figure(1); h(1).InvertHardcopy = 'off';
print('Models','-dpng');
end