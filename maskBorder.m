function A_masked = maskBorder(A,varargin)
% MASKBORDER  Mask out components on image border.
%
% Components are removed from A if a given fraction of their pixels are
% within a certain distance of the image boundary. Masked results are
% returned.
%
%   A_masked = MASKBORDER(A) masks A using default 20 pixel boundary and
%   0.9 cutoff fraction.
%
%   A_masked = MASKBORDER(A, BORDER) masks A using BORDER pixel boundary
%   and default 0.9 cutoff fraction.
%
%   A_masked = MASKBORDER(A, BORDER, FRACTION) masks A using BORDER pixel
%   boundary and FRACTION cutoff fraction.

Defaults = {20,0.9};
Defaults(1:nargin-1) = varargin;

border = Defaults{1};
border = border + 1; % for indexing
fraction = Defaults{2};
componentsToDelete = [];

for componentIdx = 1:size(A.spatial,3)
    componentImage = A.spatial(:,:,componentIdx);
    componentSize = nnz(componentImage);
    componentSizeOutsideMask = nnz(componentImage(border:end-border,border:end-border));
    componentSizeWithinMask = componentSize - componentSizeOutsideMask;
    if componentSizeWithinMask/componentSize > fraction
        componentsToDelete = [componentsToDelete; componentIdx];
    end

end

A_masked = A;
if ~isempty(componentsToDelete)
    fieldsToUpdate = ["F0";"DFF";"rawDFF";"F";"rawF"];
    for i=1:length(fieldsToUpdate)
        A_masked.(fieldsToUpdate(i))(componentsToDelete,:) = [];
    end
    A_masked.spatial(:,:,componentsToDelete)= [];
end