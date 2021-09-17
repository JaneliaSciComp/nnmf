function datafilt = whitenedMatchedFilter(data,locs, window)
N = 2*length(data)-1;

censor = false(size(data)); censor(locs(:) + window(:)') = true;
censor = censor(1:size(data,1), 1:size(data,2));

noise = data(~censor);
Nf2 = make2sided(pwelch(noise,1000,[],N))';
scaling = 1./sqrt(Nf2); %multiplier for whitening

%noiseScaled = real(ifft(fft(noise,N).*scaling));
dataScaled = real(ifft(fft(data,N).*scaling));

PTDscaled = dataScaled(locs(:)+window(:)');
PTAscaled = mean(PTDscaled,1);

datafilt = conv(dataScaled, fliplr(PTAscaled), 'same');
datafilt = datafilt(1:length(data));

    function Fout = make2sided(F)
        Fout = [F ; flipud(F(1:end-1))];
    end
end