function [fig] = plot_MVN(mvn_samples,mu,ss,free_params,i,dirName)
fig=figure(1);clf; set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 20 20], 'PaperUnits', 'centimeters', 'PaperSize', [20 20]); 
dx = .001; dy = .001; cols = get(gca,'colororder'); 

n_pars = length(mu);%length(mu(1:10)); 
count = 1; 
for kk = 1:n_pars
    for k = 1:n_pars
        subplotHJ(n_pars,n_pars,count,dy,dx); 
        hold on
        if (k>=kk | count==length(mu)*(length(mu)-1)+1)
        try
%             sc=scatter(mvn_samples(:,k),mvn_samples(:,kk),5, 'MarkerEdgeColor','none', 'MarkerFaceColor',cols(1,:),'MarkerFaceAlpha', .7);             
            
            mu0=mu([k kk]); % means
            scatter(mu0(1),mu0(2),7, '*', 'MarkerEdgeColor','k','LineWidth',.4);  
                        
            ee=error_ellipse(ss([k kk],[k kk]),mu([k kk]),'conf',.95); ee.Color = cols(1,:);    
        catch
        end
%         xlim([-5 5]); ylim([-5 5]);
        xlabel(['$log10\lambda_{', num2str([free_params(k)]), '}$'],'Interpreter','latex')
        ylabel(['$log10\lambda_{', num2str([free_params(kk)]), '}$'],'Interpreter','latex')
        if count==length(mu)*(length(mu)-1)+1
            axis off; xlim([6 8]); ylim([6 8]); 
            legend('mvnrnd samples', '\mu', 'mvnrnd 95% CI'); legend('boxoff'); legend('Location', 'northeast'); 
        end
        else
            axis off
        end        
        count = count + 1; 
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',8, 'defaultTextFontSize',8, 'FontName', 'Helvetica');

figname = [dirName,'MVNs',num2str(i)]; 
print(figname,'-depsc', '-r600'); 

end