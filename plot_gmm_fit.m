function GMmodel = plot_gmm_fit(data, nGaussians, varargin)
% PLOT_GMM_FIT Plots the Gaussian mixture model fit of data, returning
% the mixture model

plot_data = true;
if nargin==3
    plot_data = varargin{1};
end


options = statset('MaxIter',1000);
GMmodel = fitgmdist(data, nGaussians, 'options', options);
pdfPoints = (min(data):(max(data)-min(data))/1000:max(data))';

wholePDF = pdf(GMmodel,pdfPoints);
individualPDFs = zeros(length(pdfPoints), nGaussians);
for gaussianIndex = 1:nGaussians
    mu = GMmodel.mu(gaussianIndex);
    sigma = sqrt(GMmodel.Sigma(gaussianIndex));
    componentProportion = GMmodel.ComponentProportion(gaussianIndex);
    individualPDFs(:,gaussianIndex) = componentProportion*normpdf(pdfPoints, mu, sigma);
end

if plot_data
    fig = figure();
    ax = axes(fig);
    hold(ax, 'on');
    set(ax, 'FontSize', 12);

    histogram(ax, data,'Normalization','pdf','facecolor','b','facealpha',0.7);
    plot(ax, pdfPoints, wholePDF, 'r','linewidth', 4)
    plot(ax, repmat(pdfPoints, 1, nGaussians), individualPDFs, 'linewidth', 2)
    xlabel(ax, '\DeltaF/F', 'FontSize', 14);
    ylabel(ax, 'Count (PDF Normalization)', 'FontSize', 14);
end
end
