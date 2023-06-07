function A = processMiniData_edited(varargin)
addpath(pwd);

%process recordings of GluSnFR minis from the Schreiter Nikon Scope

%parameters
nIter = 5; %number of outer loops to run
bgFrac = 0.2; %fraction of extra components to use to fit non-constrained background in nmf
sparseFac = 0.15; % set to 0 all pixels with less than sparseFac of the maximum weight
maxN = 1600;
cameraOffset = 0; %(was 400)hamamatsu cameras have an offset of 100 per pixel, 400 for 2x2 binning
sigma = 2; %area to consider for sources (gaussian sigma of initial weights, full window width will be ~ 6*sigma+1)
[b1,a1] = butter(4, 0.02); %lowpass filter for defining F0/bleaching correction
[b2,a2] = butter(4, 0.08, 'high'); % 0.08 highpass filter for generating correlation image

hFvs = [];
% %load files
% if ~nargin
%     dr = uigetdir('\\dm11\podgorskilab\GluSnFR3');
% end
% %find all datasets in folder
% fns = dir([dr filesep '*.nd2']);
% fns = {fns.name};
% %discard all APO and SAT images
% keep = cellfun(@isempty, strfind(fns, 'APO'));
% fns = fns(keep);
% keep = cellfun(@isempty, strfind(fns, 'SAT'));
% fns = fns(keep);

if nargin==0
    thisDr = fileparts(which('processMiniData'));
    [fns, dr] = uigetfile([thisDr filesep '*.*'], 'multiselect', 'on');
else
    [dr,name,ext] = fileparts(varargin{1});
    fns = [name ext];
    sigma = varargin{2};
end

if ~iscell(fns)
    fns = {fns};
end
for fnum = length(fns):-1:1
    
    logfile = dr + "/log.txt";
    logfileID = fopen(logfile,'w');
    %     if exist([dr filesep fns{fnum}(1:end-4) '.mat'], 'file')
    %         disp(['Already processed:' fns{fnum} ', skipping.'])
    %         continue
    %     end
    A = struct();%initialize output, one output struct for each input file
    % This will fail if running on cluster, so try to read in tiff
    logger(logfileID,"Reading input file")
    f = bfopen([dr filesep fns{fnum}]);
    IM = double(cat(3, f{1}{:,1}))-cameraOffset;
    IM = imgaussfilt(IM,0.5);
    %imshow(bwareaopen(a>0,10))
    
    % gfp_dir = fileparts(fileparts(dr));
    % gfp_file = dir([gfp_dir '/' '*.tif']);
    % gfp_normalized=mat2gray(imread([gfp_file.folder '/' gfp_file.name]));

    %bw = bwareaopen(imbinarize(gfp,'adaptive'),9);
    %bw = imdilate(bw,strel('diamond',3));
    %for i=1:size(IM,3)
    %    IM(:,:,i) = IM(:,:,i).*gfp_normalized;
    %end

%     gfp=imread(maskchannel);
%     bw = bwareaopen(imbinarize(gfp,'adaptive'),9);
%     bw = imdilate(bw,strel('diamond',3));
%     for i=1:size(IM,3)
%         IM(:,:,i) = IM(:,:,i).*bw;%gfp_normalized;
%     end
    
    % IMds = IM;
    T2 = size(IM,3);
    T = size(IM,3);
    it=0;
    sz = [size(IM,1) size(IM,2)];
    meanIM = mean(IM,3);
    clear f
    
    %downsample in space (we should have binned 4x at acquisition time...)
    %IM = IM(1:2:end,1:2:end,:) + IM(2:2:end,1:2:end,:) + IM(1:2:end,2:2:end,:) + IM(2:2:end,2:2:end,:);
    
    %     sz = [size(IM,1) size(IM,2)];
    %     T = size(IM,3);
    %
    %     %downsample in time
    %     IMds = IM;
    %     nIter = 0; %number of time downsampling iterations
    %     for it = 1:nIter
    %         IMds = IMds(:,:,1:2*floor(end/2));
    %         IMds=IMds(:,:,1:2:end) + IMds(:,:,2:2:end);
    %     end
    %     IMds = IMds./(2.^it); %keep scale;
    %
    %     meanIM = mean(IMds,3);
    %     T2 = size(IMds,3);

    
    %compute F0
    logger(logfile,"computing F0...");

    v = 0;
    a = 0.04;
    %compute a leaky cumulative minimum
    e1 =  medfilt2(reshape(IM, [], size(IM,3)), [1 3], 'symmetric');
    logger(logfileID,'medfilt');
    for t = 2:size(e1,2)
        e1(:,t) = min(e1(:,t), e1(:,t-1) + v + a);
        v = max(0, e1(:,t) - e1(:,t-1));
    end
    clear v a
    logger(logfileID,'done');
    %use it to make a smooth F0 that obeys the data minima
    F0 = filtfilt(b1,a1,e1');
    logger(logfileID,"filtfilt");
    for ii = 1:6
        delta = min(0,e1'-F0);
        delta([1:20, end-19:end],:) = 0;
        F0 = filtfilt(b1,a1,F0+2*delta); %2*delta to accelerate convergence
    end
    logger(logfileID,"done");

    clear e1 v a
    F0 = reshape(F0', size(IM));
    logger(logfileID,"reshaped");
    
    %compute F0 and dF
    % dF = IM-F0;
    % if it==0
    %    F0 = F0ds;
    % else
    %    F0 = reshape(interp1((0.5:1:T2).*(2.^it), reshape(F0ds,[],T2)',  1:T, 'linear', 'extrap')', size(IM));
    % end
    dF = IM - F0;
    
    %ensure that F0 is always>0; We will additionally regularize when
    %calculating DFF
    F0 = F0 - min(0, min(F0,[],3));
    clear IM
    
    %save out the downsampled tiff stacks
    logger(logfileID,"Saving downsampled movies...");
    if ~exist([dr filesep fns{fnum}(1:end-4) '_DF.tif'], 'file')
        saveastiff(dF, [dr filesep fns{fnum}(1:end-4) '_DF.tif']);
        saveastiff(F0, [dr filesep fns{fnum}(1:end-4) '_F0.tif']);
        saveastiff(dF./F0, [dr filesep fns{fnum}(1:end-4) '_DFF.tif']);
        %bfsave(uint16(IMds), [dr filesep fns{fnum}(1:end-4) '_ds.tif']);
    end

    %highpass filter for correlation image
    dFhp = permute(filtfilt(b2,a2,permute(dF, [3 1 2])), [2 3 1]);
    
    %compute correlation image on downsampled recording
    logger(logfileID, "computing correlation image");

    ss = sum(dFhp.^2,3);
    vertC = sum(dFhp .* circshift(dFhp, [1 0 0]),3)./sqrt(ss.*circshift(ss, [1 0 0]));
    horzC = sum(dFhp .* circshift(dFhp, [0 1 0]),3)./sqrt(ss.*circshift(ss, [0 1 0]));
    C = nanmean(cat(3, horzC, circshift(horzC,1,2), vertC, circshift(vertC, 1,1)),3);
    
    C(isnan(C))=0; % TODO: is this ok??
    
    A.corrIM = C;
    logger(logfileID,"done");
    
    %find peaks in correlation image to initialize cNMF
    logger(logfileID,"finding peaks");

    C = imgaussfilt(C,0.5);
    C2 = C;
    Cthresh = 2.*(median(C(:))-prctile(C(:),1));
    C2(C2<(median(C(:))+Cthresh)) = 0;
    BW = imregionalmax(C2);
    C2(imdilate(BW, strel('disk', 6*sigma))) = 0;
    while any(C2(:))
        BW = BW | imregionalmax(C2);
        C2(imdilate(BW, strel('disk', 6*sigma))) = 0;
    end
    BWinds = find(BW(:));
    if length(BWinds)>maxN
        BWvals = C(BWinds);
        [~,sortorder] = sort(BWvals, 'descend');
        BWinds = BWinds(sortorder(1:maxN));
    end
    
    nComp = length(BWinds);
    [Pr,Pc] = ind2sub(sz, BWinds); %locations of putative release sites
    logger(logfileID,"done");
    
    %global cNMF to get footprints and traces for each mini location, with
    %mild unmixing of overlapping signals
    nBG = floor(nComp.*bgFrac+10);
    W0 = zeros([sz nComp+nBG]); %start with twice as many components as local maxima
    W0(sub2ind(size(W0), Pr',Pc', 1:nComp)) = 1; %initialize one pixel for every component
    W0 = imgaussfilt(W0,0.5)+imgaussfilt(W0,sigma, 'FilterSize', 6*ceil(sigma)+1);
    W0 = reshape(W0, prod(sz),nComp+nBG);
    W0(:, nComp+1:end) = rand(size(W0,1), size(W0,2)-length(Pr))./size(W0,1); %background components are initialized random
    
    logger(logfileID,string(nComp)+" nmf");

    %Use multiplicative updates NMF, which makes it easy to zero out pixels
    opts1 = statset('MaxIter', 20,  'Display', 'iter');%, 'UseParallel', true);
    [W0,H0] = nnmf(reshape(dF,[],T2), nComp+nBG,'algorithm', 'mult', 'w0', W0, 'options', opts1); %!!nnmf has been modified to allow it to take more than rank(Y) inputs
    for bigIter = 1:nIter
        logger(logfileID,"outer loop " + int2str(bigIter) + " of " + int2str(nIter));
        nW0 = sum(W0>0,1);
        
        %apply sparsity
        setZero = W0<(sparseFac.*max(W0,[],1));
        setZero(:, nW0<=9) = false; %don't shrink any more once below 9 pixels
        W0(setZero) = 0;
        
        %apply contiguous constraint
        smallComps = nW0<(prod(sz)*0.05); %sparse components, that we will apply contiguous constrains to; we have to do this because matlab's nnmf reorders components
        W0 = reshape(W0, sz(1),sz(2),[]);
        for comp = find(smallComps)
            [maxval, maxind] = max(reshape(W0(:,:,comp),1,[]));
            [rr,cc] = ind2sub(sz, maxind);
            if nW0(comp)<5
                W0(max(1,min(end,rr+(-1:1))),max(1,min(end,cc+(-1:1))),comp) = maxval/3;
            end
            W0(:,:,comp) = W0(:,:,comp).*bwselect(W0(:,:,comp)>0, cc,rr, 4);
        end
        W0 = reshape(W0, prod(sz),[]);
        [W0,H0] = nnmf(reshape(dF,[],T2), nComp+nBG,'algorithm', 'mult', 'w0', W0, 'h0', H0, 'options', opts1);
    end
    logger(logfileID,"done");
    
    %get the traces for the full-time-resolution dataset, without nonnegativity constraints
    %constraint
    dF = reshape(dF,[],T);
    F0 = reshape(F0,[],T);
    Hhf = W0\dF;
    
    %merge small components if they have high correlation
    logger(logfileID,"merge components");
    nW0 = sum(W0>0,1);
    smallComps = nW0<(prod(sz)*0.01);
    recalc = false;
    activityCorr = corr((Hhf-smoothdata(Hhf,2,'movmean',75))', 'type', 'Spearman');
    for c1 = 1:size(W0,2)
        for c2 = (c1+1):size(W0,2)
            if all(smallComps([c1 c2])) && (activityCorr(c1,c2)>0.50) && any(imdilate(reshape(W0(:,c1),sz), ones(3)) & reshape(W0(:,c2),sz),'all')
                W0(:,c1) = W0(:,c1) + W0(:,c2);
                W0(:,c2) = 0;
                activityCorr(c2,:) = 0; activityCorr(:,c2) = 0;
                recalc = true;
            end
        end
    end
    if recalc
        logger(logfileID, "Some components were merged. Recalculating factorization.");
        sel = any(W0,1);
        W0 = W0(:, sel);
        H0 = H0(sel,:);
        %run some more NMF
        [W0,~] = nnmf(reshape(dF,[],T2), sum(sel),'algorithm', 'mult', 'w0', W0, 'h0', H0, 'options', opts1);
        Hhf = W0\dF; %solve for full speed data
    end
    logger(logfileID,"done");

    %reconstruct the movie and look at residuals; could be used to refine
    %component definition
    %recon = W0*Hhf;
    
    %Select only sparse components
    nW0 = sum(W0>0,1);
    smallComps = nW0<(prod(sz)*0.01);
    W0 = W0(:,smallComps);
    Hhf = Hhf(smallComps,:);
    
    %visualize components
    delete(hFvs)
    hFvs = visualize_comps(W0,sz);
    set(hFvs, 'name', fns{fnum});
    drawnow;
    
    fprintf(logfileID,'%s\n','sort components');
    %sort components by high frequency power
    ac = sum(W0,1)'.*Hhf;
    [~,sortorder] = sort(sum((ac-smoothdata(ac,2,'movmean',75)).^2,2), 'descend');
    W0 = W0(:,sortorder);
    Hhf = Hhf(sortorder,:);
    
    %compute DFF using the most sensitive pixels within the ROI
    DFF=nan(size(Hhf)); rawDFF = DFF;
    F=nan(size(Hhf)); rawF = F; Fzero = F;
    lambda = 10; %prctile(meanIM(meanIM(:)>0),1); %regularizer, larger favors selecting brighter pixels
    for comp = 1:size(Hhf,1)
        support = find(W0(:,comp)>0);
        F(comp,:) = sum((W0(support,comp).*Hhf(comp,:)),1);
        pxDFF = (W0(support,comp).*Hhf(comp,:))./(F0(support,:)+lambda);
        [~, sortorder] = sort(sqrt(sum(pxDFF.^2,2)), 'descend');
        selpix = sortorder(1:min(end,9)); %take the 9 highest-DFF pixels
        Fzero(comp,:) = sum(F0(support(selpix),:),1);
        DFF(comp,:)= sum(W0(support(selpix),comp)*Hhf(comp,:),1)./Fzero(comp,:);
        rawF(comp,:) = sum(dF(support(selpix),:),1);
        rawDFF(comp,:) = sum(dF(support(selpix),:),1)./Fzero(comp,:);
    end
    logger(logfileID,"done");

    
    A.fn = [dr filesep fns{fnum}];
    A.lambda = lambda;
    A.IM = meanIM;
    A.F0 = Fzero;
    A.DFF = DFF; A.rawDFF = rawDFF;
    A.F = F; A.rawF = rawF;
    A.spatial = reshape(W0, [sz size(W0,2)]);
    A.dr = dr;
    A.filename = fns{fnum};
    
    create_clustergram(DFF);
    create_pca_analysis_figures(DFF);
    create_dimensionality_reduction_colored_movie(A);
    save([dr filesep fns{fnum}(1:end-4) 'v2'], 'A', '-v7.3');
    clear F0 dF A
end
end

function logger(logfileID, msg)
    [tmp pid] = system('pgrep MATLAB');
    [tmp mem_usage] = system(['cat /proc/' strtrim(pid) '/status | grep VmSize']);
    mem_usage = sprintf("%i MB\n", round(str2num(strtrim(extractAfter(extractBefore(mem_usage, ' kB'), ':'))) / 1000));
    msg = string(datetime("now")) + " " + msg + ": " +  mem_usage;
    display(msg)
    fprintf(logfileID,'%s',msg);

end


function hF = visualize_comps(S, sz)
nS = size(S,2);
RGB = rand(3,nS).^2;
RGB = RGB./repmat(sum(RGB,1), 3,1);
S_RGB = sqrt([S*RGB(1,:)' S*RGB(2,:)' S*RGB(3,:)']);
S_RGB = 1.5* S_RGB./max(S_RGB(:));
S_RGB = reshape(full(S_RGB), [sz 3]);
hF = figure('Name', 'NMF components'); imshow(S_RGB);
end

