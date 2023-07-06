function summarizeMiniData
    
thisDr = fileparts(which('summarizeMiniData'));
[fns, dr] = uigetfile([thisDr filesep '*.mat'], 'multiselect', 'on');
if ~iscell(fns)
    fns = {fns};
end

    %we'll save data in this folder:
    savedr = [dr filesep 'SummaryData'];
    if ~exist(savedr,'dir')
        mkdir(savedr);
    end
    
    
    %we will output:
    B.consNum = nan(1,length(fns)); 
    B.wellNum = nan(1, length(fns));
    B.fovNum = nan(1, length(fns));
    
    B.fileName = cell(1,length(fns));
    B.traces = cell(1,length(fns)); %spikes
    %B.S = cell(1,length(fns)); %spikes
    %B.C = cell(1,length(fns)); %deconvolved trace
    B.tau = cell(1,length(fns)); %time constant
    B.snr = cell(1,length(fns)); %SNR
    B.freq = cell(1,length(fns)); % Event rate
    
    for fnum = length(fns):-1:1
        disp(['Processing: ' int2str(fnum)])
%         pat = "C" + digitsPattern;
%         Cmatch = extract(fns{fnum}, pat);
%         B.consNum(fnum) = str2double(Cmatch{1}(2:end));
%         pat = "W" + digitsPattern(1);
%         Wmatch = extract(fns{fnum}, pat);
%         B.wellNum(fnum) = str2double(Wmatch{1}(2));
%         pat = "F" + digitsPattern;
%         Fmatch = extract(fns{fnum}, pat);
%         B.fovNum(fnum) = str2double(Fmatch{1}(2));
        
        load([dr filesep fns{fnum}], 'A');
        [S,C, tau, snr, freq, inds, traces] = runFOOPSI(A.DFF');
        
        B.fileName{fnum} = fns{fnum};
        
        %we should probably filter out garbage traces at this point??
        %B.S{fnum} = S;
        %B.C{fnum} = C;
        B.traces{fnum} = single(traces);
        B.inds{fnum} = inds;
        B.tau{fnum} = tau;
        B.snr{fnum} = snr;
        B.freq{fnum} = freq;
        
        %summary statistics for this recording
        selFast = tau<18;
        selSNR = snr>3;
        selFreq = freq>(5/6000);
        B.medTau(fnum) = median(tau(selSNR & selFreq));
        B.medSNR(fnum) = median(snr(selFast & selFreq));
        B.medFreq(fnum) = median(freq(selFast & selSNR));
    end
    
    save([savedr filesep 'minidata'], 'B', '-v7.3');
end

function [S, C, tau, snr, freq, inds, traces] = runFOOPSI(traces)
%use OASIS (Friedrich et al, PLoS Comp. Bio. 2017) to estimate AR parameters and detect events
%we use unconstrained deconvolution (no sparsity penalty) to be able to judge the
%separation of signal from noise, and report d_prime

options.type = 'ar1';
options.optimize_b = false;
options.optimize_pars = false; %true
options.method = 'foopsi'; %'foopsi', 'constrained', 'thresholded', 'mcmc'
%options.max_tau = 9; %for ar(1) only
options.tau_range = [9 9]; %for ar(1) only

S = nan(size(traces)); %spike matrix
C =  nan(size(traces)); %deconvolved signal
inds = false(size(traces));
tau = nan(1, size(traces,2)); %decay constant
snr = nan(1, size(traces,2)); %snr
freq = nan(1, size(traces,2)); %event frequency (per frame)
[b2,a2] = butter(4, [0.01 0.5], 'bandpass');


%compute F0
[b1,a1] = butter(4, 0.004);
v = 0;
a = 1e-5;
nIter = 10;
%compute a leaky cumulative minimum
e1 =  medfilt2(traces, [7 1], 'symmetric');
for t = 2:size(traces,1)
    e1(t,:) = min(e1(t,:), e1(t-1,:) + v + a);
    v = max(0, e1(t,:) - e1(t-1,:));
end
%use it to make a smooth F0 that obeys the data minima
F0 = filtfilt(b1,a1,e1);
for iter = 1:nIter 
    delta = min(0,e1-F0);
    delta([1:51, end-51:end],:) = 0;
    F0 = filtfilt(b1,a1,F0+2*delta);
end

for traceIx = 1:size(traces,2)
    trace = traces(:,traceIx)-F0(:,traceIx);

   [c, s, optsOut] = deconvolveCa(trace, options);
   S(:,traceIx) = s;
   C(:,traceIx) = c;
   
   
   %filter for event detection
   HP = filtfilt(b2,a2,c);
   traceStd = sqrt(estimatenoise(trace));
   [pks, locs] = findpeaks(HP);
   
   sortedPks  = sort(pks, 'descend');
   if length(sortedPks)<5
      snr(traceIx) = 0;
   else
      snr(traceIx) = max(0, sortedPks(5)./traceStd);
   end
   
   %compute frequency
   sel = pks>3*traceStd;
   locs = locs(sel);
   inds(locs,traceIx) = 1; 
   freq(traceIx) = length(locs)./length(trace);
   
   %snr(traceIx) = max(0, (var(c)./estimatenoise(trace))-1);
   tau(traceIx) = -1/log(optsOut.pars);
   
   %nonmax suppression/event detection
   %Sfilt = filtfilt(b2,a2,c);
%    Sfilt = smoothdata(s, 'movmean',2); 
%    [x, locs] = findpeaks(Sfilt);
%    thresh = 3*median(x);
%    sel = x>thresh;
%    inds(locs(sel),traceIx) = 1; 
%    freq(traceIx) = sum(sel)./length(trace);

   traces(:,traceIx) = trace;
end
end






