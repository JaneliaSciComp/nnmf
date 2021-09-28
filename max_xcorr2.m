function D2 = max_xcorr2(ZI, ZJ)
%MAX_XCORR2 Caclulates the maximum cross correlation

D2 = zeros(size(ZJ,1),1);
for i=1:size(ZJ,1)
    [corr, lags] = xcorr(ZI, ZJ(i,:));
    [max_corr, max_corr_idx] = max(corr);
    D2(i) = abs(lags(max_corr_idx));
    %D2(i)=corr(max_corr_idx);
end
%[cross_correlation, lags] = xcorr2(ZI,ZJ);
%[~, D2] = max(cross_correlation, [], 2);