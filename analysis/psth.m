function psth(freqs,binsize,bins,tLength,nTrials,fh)
% spikes: vector of spike times in seconds
% binsize: Bin size in seconds
% tLength: Total recorded time
% nTrials: The number of trials (or cells) to average over
% fh: Figure handle

% bins = binsize/2 : binsize : (tLength - binsize/2); % specify centres for hist(), or edges for histc()
% freqs = hist(spikes,bins)/(nTrials*binsize); % Normalise
%freqs = freqs - mean(freqs);
freqs = freqs/(nTrials*binsize); % Normalise
figure(fh); bar(bins*1000,freqs,1); % Plot as msec
xlim([0,tLength*1000]); %+(step/2)
xlabel('Time (ms)'); %Seconds (Binned)
%ylabel('Average Spikes/s'); % Without normalising would just be frequency
ylabel('Spikes/s');

% Grun & Rotter 2010 Ch.2
% http://www.ton.scphys.kyoto-u.ac.jp/~shino/toolbox/english.htm

% optN = sshist(x); hist(x,optN);
% optW = sskernel(x); ksdensity(x,'width',optW);