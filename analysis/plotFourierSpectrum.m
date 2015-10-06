function plotFourierSpectrum(freqs,binsize,fh)%tLength,fh)
% spikes: vector of spike times in seconds
% binsize: Bin size in seconds
% tLength: Total recorded time
% fh: Figure handle

% http://www.mathworks.com/products/matlab/demos.html?file=/products/demos/shipping/matlab/sunspots.html

% bins = binsize/2 : binsize : (tLength - binsize/2); % specify centres for hist(), or edges for histc()
% freqs = hist(spikes,bins);

%freqs = freqs - mean(freqs);

% bins = 0:binsize:tLength;
% freqs = histc(spikes,bins);
% freqs(end) = [];

Fs = 1/binsize; % Sampling frequency (Hz)
NyquistLim = Fs/2; % (1/sampInt)/2
disp(['Nyquist Limit = ',num2str(NyquistLim,3),' Hz']);

NFFT = 2^nextpow2(length(freqs)); % Next power of 2 from length of freqs %http://blinkdagger.com/matlab/matlab-fft-and-zero-padding/
FFT = fft(freqs,NFFT);%/length(freqs);%Fs;
FFT(1) = []; % This is simply the sum of the data and can be removed
% Plot first half of FFT (second half is a reflection of the first)
f = NyquistLim * linspace(0,1,NFFT/2 + 1); % Gives f up to Nyquist limit NyquistLim * (0:NFFT-1); %
%figure(fh); plot(f,2*abs(FFT(1:NFFT/2 +1))/NFFT);    ylabel('|Amplitute|');
figure(fh); plot(f,(abs(FFT(1:NFFT/2 + 1)).^2)/NFFT);   ylabel('Power'); % Periodogram
xlabel('Frequency (Hz)');
xlim([0,NyquistLim]); %xlim([-1,NyquistLim+1]);
