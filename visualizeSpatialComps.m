function S_RGB = visualizeSpatialComps(A)
% VISUALIZESPATIALCOMPS  Create an interactive visualization of spatial
% components of A, allowing clicking on the components to display their
% fluorescence traces. Returns the spatial component colored image.

S = A.spatial;
[maximumIntensity,maximumIntensityCluster] = max(S,[],3);

function ClearLinesFromAxes()
  axesHandlesToChildObjects = findobj(gca, 'Type', 'line');
  if ~isempty(axesHandlesToChildObjects)
    delete(axesHandlesToChildObjects);
  end  
  return; % from ClearLinesFromAxes
end

function mytestcallback(~,~)
    pt = get(gca,'CurrentPoint');
    x = floor(pt(1,1));
    y = floor(pt(1,2));

    if maximumIntensity(y,x)>0
        ClearLinesFromAxes()
        clusterIdx = maximumIntensityCluster(y,x);
        subplot(1,2,1)
        [B,~] = bwboundaries(A.spatial(:,:,clusterIdx)>0);
        for k = 1:length(B)
           boundary = B{k};
           plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
        end
        subplot(1,2,2)
        plot(A.DFF(clusterIdx,:))  
        title("Cluster " + clusterIdx)
        xlabel('Frame')
        ylabel('{\Delta}F/F')
    end
end

sz = size(S);
nS = sz(3);
S = reshape(S, [], nS);
RGB = rand(3,nS).^2; RGB = RGB./repmat(sum(RGB,1), 3,1);
S_RGB = sqrt([S*RGB(1,:)' S*RGB(2,:)' S*RGB(3,:)']);
S_RGB = 1.5* S_RGB./max(S_RGB(:));
S_RGB = reshape(full(S_RGB), [sz(1:2) 3]);
fig = figure();
subplot(1,2,1);
imshow(S_RGB)
hold on;
subplot(1,2,2)
xlabel('Frame')
ylabel('{\Delta}F/F')
xlim([0 size(A.DFF,2)])
set(fig,'WindowButtonDownFcn',@mytestcallback)
end
