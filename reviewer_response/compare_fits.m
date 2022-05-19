%% load data and fix rng for reproducibility
close all; clear all;
data = load('peaks_data.mat').peaks11_compiled';

%% try fitting the data using 1-4 gaussians and pick best based on lowest AIC (https://en.wikipedia.org/wiki/Akaike_information_criterion)
gmm = cell(4);
AIC = zeros(4,1);
options = statset('MaxIter',1000);
num_trials = 1000;
track_best_num_gaussians = zeros(num_trials,1);
% try over a range of trials
rng(1)
for trial = 1:num_trials
    for k = 1:4
        gmm{k} = fitgmdist(data, k, 'options', options);
        AIC(k)= gmm{k}.AIC;
    end
    [best_aic, best_num_gaussians] = min(AIC);
    % probability_density_gmm = gmm{best_num_gaussians};
    track_best_num_gaussians(trial) = best_num_gaussians;
end
best_num_gaussians = mode(track_best_num_gaussians);
fprintf('Best number of guassians over %d trials: %d.\n\n', num_trials, best_num_gaussians);
%% compare best gmm to lognormal
probability_density_lognormal = fitdist(data,'Lognormal');

% plot pdfs
track_best_model = zeros(1,num_trials);
for trial=1:num_trials
    probability_density_gmm = plot_gmm_fit(data,best_num_gaussians,trial==num_trials);
    pdf_points = (min(data):(max(data)-min(data))/1000:max(data))';
    y_lognormal = pdf(probability_density_lognormal, pdf_points);

    if trial==num_trials
        plot(pdf_points, y_lognormal,'k','linewidth', 4)
        legend('Empirical Distribution','GMM Fit','Gaussian 1', 'Gaussian 2', 'Gaussian3','Lognormal Fit','location','best' )
    end

    % calculate AICs
    loglikelihood_gmm = -probability_density_gmm.NegativeLogLikelihood;
    loglikelihood_lognormal = -lognlike(probability_density_lognormal.ParameterValues,data);

    AIC_gmm = aicbic(loglikelihood_gmm,8); % equivalent to probability_density_gmm.AIC (https://stats.stackexchange.com/questions/229293/the-number-of-parameters-in-gaussian-mixture-model)
    AIC_lognormal = aicbic(loglikelihood_lognormal, 2);
    track_best_model(trial) = AIC_gmm<AIC_lognormal;
end

model_names = {'Lognormal'; 'GMM'};
fprintf('Best model over %d trials is %s.\n\n', num_trials, model_names{mode(track_best_model)+1});

% Print loglikelihood
fprintf('Loglikelihood for GMM: %e.\n', loglikelihood_gmm);
fprintf('Loglikelihood for Lognormal: %e.\n\n', loglikelihood_lognormal);

% Print AICs
fprintf('AIC for GMM: %e.\n', AIC_gmm);
fprintf('AIC for Lognormal: %e.\n\n', AIC_lognormal);

% plot cdfs
fig = figure();
ax = axes(fig);
hold(ax, 'on');
set(ax, 'FontSize', 12);
[f,x] = ecdf(data);
plot(x,f,'b','linewidth',3)
plot(pdf_points, probability_density_gmm.cdf(pdf_points), 'g','linewidth',1.5);
plot(pdf_points, probability_density_lognormal.cdf(pdf_points),'r','linewidth',1.5);
xlabel('\DeltaF/F', 'FontSize', 14);
ylabel('CDF', 'FontSize', 14);
legend('Empirical CDF', 'GMM CDF', 'Lognormal CDF','location','best' )

% plot qq-plots

% taken from qqplot.m
y = probability_density_lognormal;
[sorted_x,~] = sort(data);
probabilities = plotpos(sorted_x);
icdf_xs_lognormal = icdf(probability_density_lognormal, probabilities);
icdf_xs_gmm = gmm_icdf(probability_density_gmm, probabilities, data);

fig = figure();
ax = axes(fig);
hold(ax, 'on');
set(ax, 'FontSize', 12);
plot([min(data),max(data)],[min(data), max(data)],'b','linewidth',3);
plot(sorted_x, icdf_xs_gmm,'r+')
plot(sorted_x, icdf_xs_lognormal,'g+')
ylabel('Empirical Quantiles', 'FontSize', 14);
xlabel('Model Quantiles', 'FontSize', 14);
legend('y=x','GMM','Lognormal','location','best')
gmm_rms = mean((icdf_xs_gmm-sorted_x).^2);
lognormal_rms = mean((icdf_xs_lognormal-sorted_x).^2);

% Print Q-Q plot r-squared
fprintf('RMS quantile devation of GMM from y=x: %e.\n',gmm_rms);
fprintf('RMS quantile devation of Lognormal from y=x: %e.\n\n',lognormal_rms);
%% taken from qqplot.m
function pp = plotpos(sx)
%PLOTPOS Compute plotting positions for a probability plot
%   PP = PLOTPOS(SX) compute the plotting positions for a probability
%   plot of the columns of SX (or for SX itself if it is a vector).
%   SX must be sorted before being passed into PLOTPOS.  The ith
%   value of SX has plotting position (i-0.5)/n, where n is
%   the number of rows of SX.  NaN values are removed before
%   computing the plotting positions.

[n, m] = size(sx);
if n == 1
    sx = sx';
    n = m;
    m = 1;
end

nvec = sum(~isnan(sx));
pp = repmat((1:n)', 1, m);
pp = (pp-.5) ./ repmat(nvec, n, 1);
pp(isnan(sx)) = NaN;
end

%% get cdf estimate from gmm
function icdf_estimate = gmm_icdf(pd_gmm, p, data)
x_temp = (-max(data):1E-7:max(data))';
cdf_estimate = pd_gmm.cdf(x_temp);
icdf_estimate = zeros(size(p));
for p_index = 1:length(p)
    [~,i] = min(abs(cdf_estimate-p(p_index)));
    icdf_estimate(p_index) = x_temp(i);
end
end