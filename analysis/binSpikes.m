function [freqs,bins] = binSpikes(spikes, binsize, tLength)
%spikes(spikes(:)==0)=[];
spikes=spikes(:)'; % Converts matrix to vector % reshape(A, 1, prod(size(A)))
%spikes(isnan(spikes(:)))=[];
spikes(spikes(:)<0)=[]; % Not strictly necessary since bins start at 0
bins = 0:binsize:tLength;
freqs = histc(spikes,bins);
freqs(end) = [];
bins(end) = [];
%freqs = freqs/(nTrials*binsize); % Normalise ???