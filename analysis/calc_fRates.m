function fRates = calc_fRates(spikeTimes, nStimuli, nTransPS, nCells, transTS)
fRates = zeros(nTransPS, nCells, nStimuli);
edges = (0 : (nTransPS * nStimuli)) * transTS; % Calculate bin edges
if isvector(spikeTimes) % Pad into matrix
    spikeTimes = [spikeTimes, (ones(length(spikeTimes),1)*-1)]; % This forces histc to bin spikes for each cell even when there is at most one spike per cell
end
%spikeTimes(spikeTimes==0) = -1; % Don't count zeros padding the matrix
% No longer necessary since spikes are padded with -1
freqs = histc(spikeTimes',edges); %n(k) counts the value x(i) if edges(k) <= x(i) < edges(k+1). 
%The last bin counts any values of x that match edges(end).
%freqs(end-1,:)=freqs(end-1,:)+freqs(end,:); % include spikes at last timestep - not necessary since using C numbering last possible spike is edges(end)-1
freqs(end,:) = []; % remove n^th +1 bin (spikeTimes==edges(end)) %freqs = freqs(1:end-1,:); 
%freqs = freqs * 1/transP_Test; %-timewin makes this redundant
for s=1:nStimuli
    fRates(:,:,s) = freqs(1+(s-1)*nTransPS:s*nTransPS,:);
end