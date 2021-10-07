function [pks, locs, widths, proms] = get_dff_peaks(DFF, varargin)
%GET_DFF_PEAKS Obtains peaks from DFF trace given an optinal
%minPeakProminence, with default = 0.03.


if isempty(varargin)
    minPeakProminence=0.03;
else
    minPeakProminence = varargin{1};
end

pks = cell(size(DFF,1),1);
locs = cell(size(DFF,1),1);
widths = cell(size(DFF,1),1);
proms = cell(size(DFF,1),1);
imageData = zeros(size(DFF));

for i=1:size(DFF,1)
    [pks{i},locs{i},widths{i},proms{i}] = findpeaks(DFF(i,:),'MinPeakProminence',minPeakProminence);
    imageData(i,locs{i}) = 1;
end

imshow(imageData);
xlabel("Frame");
ylabel("Cluster");