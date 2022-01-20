function plot_cross_validation_results_1cov(CV)

set(0,'defaultTextInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')

if CV.setup.penshape == 0
    figure
    hold on; box on
    plot(CV.results.lambda,CV.results.predictive_likelihood,'-o','color',[1 1 1]*0.6)
    plot(CV.results.lambda,mean(CV.results.predictive_likelihood,1),'k-','LineWidth',2)
    set(gca,'xscale','log')
    xl=xlim;
    plot(xl,0*xl+CV.results.predictive_likelihood_accepted,'r--','LineWidth',2)
    scatter(CV.results.lambda_optimal,CV.results.predictive_likelihood_min,'ro','filled')
    xlabel('Roughness coefficient $\lambda_{\sigma}$')
    ylabel('Negative log predictive likelihood')
    zoom on
else
    
    PLmean=squeeze(mean(CV.results.predictive_likelihood,1));
    PLrange=squeeze(max(CV.results.predictive_likelihood,[],1)-min(CV.results.predictive_likelihood,[],1));
    Lmin=min(log10(CV.results.lambda));
    Lmax=max(log10(CV.results.lambda));
    dL=log10(CV.results.lambda(2))-log10(CV.results.lambda(1));

    figure
    s=surf(log10(CV.results.lambda),log10(CV.results.lambda),PLmean);
    s.FaceColor='interp';
    hold on
    scatter3(log10(CV.results.lambda_optimal(1)),log10(CV.results.lambda_optimal(2)),CV.results.predictive_likelihood_min,'ro','filled')
    xlim([Lmin Lmax])
    ylim([Lmin Lmax])
    ax1=gca;
    ax1.YTick=ax1.XTick;
    ticklab=cell(length(ax1.XTick),1);
    for i=1:length(ax1.XTick)
        ticklab{i}=['$10^{' num2str(ax1.XTick(i)) '}$'];
    end
    ax1.XTickLabel=ticklab;
    ax1.YTickLabel=ticklab;
    xlabel('Roughness coefficient $\lambda_{\sigma}$')
    ylabel('Roughness coefficient $\lambda_{\xi}$')
    zlabel('Mean predictive negative log likelihood')
    
    figure
    nanimage(log10(CV.results.lambda),log10(CV.results.lambda),PLrange)
    title('Range: Negative log predictive likelihood')
    C=colorbar;
    C.TickLabelInterpreter='latex';
    axis([Lmin-dL/2 Lmax+dL/2 Lmin-dL/2 Lmax-dL/2])
    ax1=gca;
    ax1.YTick=ax1.XTick;
    ticklab=cell(length(ax1.XTick),1);
    for i=1:length(ax1.XTick)
        ticklab{i}=['$10^{' num2str(ax1.XTick(i)) '}$'];
    end
    ax1.XTickLabel=ticklab;
    ax1.YTickLabel=ticklab;
    xlabel('Roughness coefficient $\lambda_{\sigma}$')
    ylabel('Roughness coefficient $\lambda_{\xi}$')
end

end