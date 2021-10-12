function GMmodel = plot_gmm_fit(data, nGaussians)
    % PLOT_GMM_FIT Plots the Gaussian mixture model fit of data, returning
    % the mixture model
    
    fig = figure();
    ax = axes(fig);
    hold(ax);
    set(ax, 'FontSize', 12);

    histogram(ax, data,'Normalization','pdf');

    GMmodel = fitgmdist(data, nGaussians);
    pdfPoints = (min(data):(max(data)-min(data))/1000:max(data))';
    
    wholePDF = pdf(GMmodel,pdfPoints);
    individualPDFs = zeros(length(pdfPoints), nGaussians);
    for gaussianIndex = 1:nGaussians 
        mu = GMmodel.mu(gaussianIndex);
        sigma = sqrt(GMmodel.Sigma(gaussianIndex));
        componentProportion = GMmodel.ComponentProportion(gaussianIndex);
        individualPDFs(:,gaussianIndex) = componentProportion*normpdf(pdfPoints, mu, sigma);
    end
    
    plot(ax, pdfPoints, wholePDF, 'linewidth', 4)
    plot(ax, repmat(pdfPoints, 1, nGaussians), individualPDFs, 'linewidth', 2)
    xlabel(ax, 'Bin', 'FontSize', 14);
    ylabel(ax, 'Count (PDF Normalization)', 'FontSize', 14);
end
