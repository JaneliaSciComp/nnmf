function peaks = create_dff_derivative_plot(DFF)
%CREATE_DFF_DERIVATIVE_PLOT Plot DFF derivatives
peaks = cell(size(DFF,1),1);

% TODO: Use central difference
derivative = diff(DFF')';
for i = 1:size(derivative,1)
    peaks{i} = findpeaks(derivative(i,:),...,
        'MinPeakWidth',4,...
        'MinPeakDistance',10,...
        'MinPeakProminence',0.01);
end

end

