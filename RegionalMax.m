function RegionalMax

%parameters

baseline = 0; %camera baseline
sigma =  3;

%user selects parent folder for batch processing
dr = uigetdir('C:\Users\podgorskik\Desktop\tmp\AB');

%find all tif files in folder
fns = dir([dr filesep '*.tif']);
fns = {fns.name};

savedr = [dr filesep 'outputData'];
if ~isfolder(savedr)
    mkdir(savedr);
end

%for every file...
for ii = 1:length(fns)
    %load the tiff file
    disp(['loading file:' fns{ii}])
    D = tiffread2([dr filesep fns{ii}]);
    IM = double(cell2mat(reshape({D.data}, 1,1,[])))-baseline;
    
    %despeckle
    IMd = medfilt3(IM, [3 3 1]);
    clear IM
    
    %gaussian blur
    IMd = imgaussfilt(IMd, sigma*[1 1]);
    
    BW = imregionalmax(IMd,8);
    saveastiff(uint8(BW),[savedr '/test.tif']);
end

end

