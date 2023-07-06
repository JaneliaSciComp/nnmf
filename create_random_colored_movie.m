function create_random_colored_movie(results_struct)
%CREATE_MOVIE_USING_PCA_COLORS Create movie where components are colored
%based on PCA components and intensity is based on DFF.
spatial = results_struct.spatial;
DFF = results_struct.DFF;

S=reshape(spatial,[size(spatial,1)*size(spatial,2),size(spatial,3)]);
S=S>0;
nS = size(spatial,3);
sz=[size(spatial,1),size(spatial,2)];

rescaled = mat2gray(DFF,[0 .25]);


[~,filename,~] = fileparts(results_struct.filename);
RGB = rand(3,nS).^2; 
RGB = RGB./repmat(sum(RGB,1), 3,1);

output_name = [results_struct.dr filesep filename 'RandomRecolored.tiff'];
for frame=1:size(DFF,2)

    S_scaled=S.* rescaled(:,frame)';
    S_scaled_RGB = 255*[S_scaled*RGB(1,:)' S_scaled*RGB(2,:)' S_scaled*RGB(3,:)']/255;%sqrt([S_scaled*RGB(1,:)' S_scaled*RGB(2,:)' S_scaled*RGB(3,:)']);
    S_scaled_RGB = 255*reshape(full(S_scaled_RGB), [sz 3]);
    S_scaled_RGB = uint8(S_scaled_RGB);
    if frame==1
        imwrite(S_scaled_RGB, output_name)
    else
        imwrite(S_scaled_RGB, output_name, 'writemode', 'append')
    end
end

