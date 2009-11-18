function SpikeAnalysis(rdir,simrun)

if nargin == 1
    %rdir = '~/Results/SpikeNet/';
    content = dir(rdir);
    for d = 1:length(content)
        simrun = content(d).name;
        if (content(d).isdir)
            ignoreDir = any(strcmpi(simrun, {'private','CVS','.','..'}));
            ignorePrefix = any(strncmp(simrun, {'@','.'}, 1));
            if ~(ignoreDir || ignorePrefix)
                SpikeAnalysis(rdir,simrun);
            end
        end
    end
    return;
end

%path(path, '~/Code/Spike/');
eval(['cd ',strcat(rdir,simrun)]);
gspath='/usr/local/bin/gs';

%%
defaults;
parameters;
%eval(['! cp ',wdir,'*.dat ',mdir,'.']);
%%
close all;
% Load presynaptic connections for excitatory cells
%PRESYNCNX_EE = import_presyncnx('presyncnx_EE.dat') + 1; % Since C starts at 0

% if isempty(rfpath)
%     pfile = 'outputs';
% else
%     [token, remain] = strtok(rfpath,'/');
%     pfile = strtok(remain(2:end),'.');
% end

%% Load all spike time data

plotInhib=0;
outLayer = nLayers - 1;
% Place into cells?
if pretrain == 1
    ptESpikes = cell(1,nLayers);
end
ESpikes = cell(1,nLayers);
ISpikes = cell(1,nLayers);
for l=0:outLayer
    prefix = strcat(['L',int2str(l)]);
    eval(['load ',prefix,'affNeuronsEI.dat']);
    eval(['load ',prefix,'affNeuronsElE.dat']);
    eval(['load ',prefix,'affNeuronsIE.dat']);
    eval(['load ',prefix,'affNeuronsII.dat']);
    if (l>0)
        eval(['load ',prefix,'affNeuronsEfE.dat']);
        eval(['load ',prefix,'weightsEfE.dat']);
    end
    if pretrain == 1
        ptESpikes{l+1} = dlmread(strcat('pt',prefix,'ExcitSpikes.dat'));
        ptESpikes{l+1}(:,1) = []; % Delete first column
    end
    ESpikes{l+1} = dlmread(strcat(prefix,'ExcitSpikes.dat'));
    ESpikes{l+1}(:,1) = []; % Delete first column
    ISpikes{l+1} = dlmread(strcat(prefix,'InhibSpikes.dat'));
    ISpikes{l+1}(:,1) = []; % Delete first column
end



%% Information theory analysis

% bpath = '~/bin/analysis/Release/';
sLayer = nLayers-1; % Selected layer equal to output layer in 'C' numbering
% maxBins = nTransPS;
% minTrans = 4; %nTransPS;
% maxTrans = nTransPS;
% decoding = 1;
nBestCells = 5;

% Bin the Excitatory output spikes
if pretrain == 1
    ptoutFrates = calc_fRates(ptESpikes{nLayers}, nStimuli, nTransPS, nExcit, transP_Test * 1/DT);
    ptINFO = infoAnalysis(ptoutFrates, nStimuli, nTransPS, nExcit, nBestCells, sLayer, transP_Test, simrun);
end
outFrates = calc_fRates(ESpikes{nLayers}, nStimuli, nTransPS, nExcit, transP_Test * 1/DT);
INFO = infoAnalysis(outFrates, nStimuli, nTransPS, nExcit, nBestCells, sLayer, transP_Test, simrun);
 

%% Plot results of information analysis

% Mutual single cell information
% plot(sort(msci{4},'descend')); %Sort 4th column (corrected info) and plot
figure();
plot(sort(INFO{1},'descend')); %Sort 4th column (corrected info) and plot
if pretrain == 1
    hold on
    plot(sort(ptINFO{1},'descend'),'--r');
    legend('Trained','Untrained');
    hold off
end
axis([1 nExcit 0 log2(nStimuli)]);
xlabel('Cell rank');
ylabel('Information (bits)');
title('Mutual single cell information');
saveas(gcf,'Mutual single cell information','epsc');
eps2pdf('Mutual single cell information.eps',gspath);
saveas(gcf,'Mutual single cell information','png');

% Stimulus specific single cell information
% plot(sort(sssci{2},'descend')); %Sort 2nd column and plot
figure();
plot(sort(INFO{2},'descend')); %Sort 2nd column and plot
if pretrain == 1
    hold on
    plot(sort(ptINFO{2},'descend'),'--r');
    legend('Trained','Untrained');
    hold off
end
axis([1 nExcit 0 log2(nStimuli)]);
xlabel('Cell rank');
ylabel('Information (bits)');
title('Stimulus specific single cell information');
saveas(gcf,'Stimulus specific single cell information','epsc');
eps2pdf('Stimulus specific single cell information.eps',gspath);
saveas(gcf,'Stimulus specific single cell information','png');

% Multiple cell information
% plot(mci{3}); %Plot 3rd column (already sorted ascending)
figure();
plot(INFO{3}); %Plot 3rd column (already sorted ascending)
if pretrain == 1
    hold on
    plot(ptINFO{3},'--r');
    legend('Trained','Untrained');
    hold off
end
axis([1 (nBestCells*nStimuli) 0 log2(nStimuli)]);
xlabel('Number of best cells');
ylabel('Information (bits)');
title('Multiple cell information');
saveas(gcf,'Multiple cell information','epsc');
eps2pdf('Multiple cell information.eps',gspath);
saveas(gcf,'Multiple cell information','png');

%% Perform network analysis
% Calculate weight bins and centres
wbins = 0.025:0.05:0.975;
% bstart=TEST_PERIOD_TS/(NSTIMULI*2);
% bstep=bstart*2;
% bend=TEST_PERIOD_TS-bstart;
% 
% centres=(bstart:bstep:bend);
% 
% bins=zeros(NSTIMULI,NEXCIT);
% for c=1:size(output_spikes,2)
%     tmp=output_spikes(:,c);
%     tmp(tmp(:)==0)=[];
%     bins(:,c)=hist(tmp,centres);
% end

testPeriodMS = transP_Test * nStimuli * nTransPS * 1000;
trainPeriodMS = transP_Train * nStimuli * nTransPS * 1000;

step = 0.005; % seconds
tstart = step/2;
tend = (testPeriodMS/1000)-step/2;
testBins = tstart:step:tend;

spikes = ESpikes{nLayers}.*DT; %eval(['L',int2str(outLayer),'ExcitSpikes']).*DT;
%spikes = L1ExcitSpikes.*DT;
spikes(spikes(:)==0)=[];

Fs = 1/step;

figure()
freqs = hist(spikes,testBins)/(nExcit*step); % Normalise
%freqs = freqs - mean(freqs);
bar(testBins,freqs);
xlim([0,testPeriodMS/1000]); % Redo
xlabel('Seconds (Binned)');
ylabel('Frequency');
title('Output Layer Excitatory Spike Frequency Histogram');
saveas(gcf,strcat('L',int2str(outLayer),'ESpikeHist',simrun),'epsc');
eps2pdf(strcat('L',int2str(outLayer),'ESpikeHist',simrun,'.eps'),gspath);
saveas(gcf,strcat('L',int2str(outLayer),'ESpikeHist',simrun),'png');

figure()
%freqs=freqs(1:length(freqs)/2);
freqs = freqs - mean(freqs);
NFFT = 2^nextpow2(length(freqs)); % Next power of 2 from length of freqs
FFT = fft(freqs,NFFT)/Fs;
f = linspace(0,Fs/2,NFFT/2);
plot(f,2*abs(FFT(1:NFFT/2)));
xlabel('Frequency (Hz)');
ylabel('|Amplitute|');
title('Output Layer Fourier Spectrum');
xlim([-1,(Fs/2+1)]);
saveas(gcf,strcat('L',int2str(outLayer),'ESpikeFourier',simrun),'epsc');
eps2pdf(strcat('L',int2str(outLayer),'ESpikeFourier',simrun,'.eps'),gspath);
saveas(gcf,strcat('L',int2str(outLayer),'ESpikeFourier',simrun),'png');

%% Plot Rasters and weights

% Plot pre-training results
if pretrain==1
    figure();
    pt_output_spikes = dlmread('ptL1ExcitSpikes.dat');
    load ptL1weightsEfE.dat;
    subplot(2,3,1);
    %imagesc(sort_matrix(ptL1weightsEfE'));
    imagesc(ptL1weightsEfE');
    xlabel('Input neuron');
    ylabel('Output neuron'); % Check these labels
    title('Synaptic Weight Matrix');

    subplot(2,3,4);
    hist(reshape(ptL1weightsEfE,1,[]),wbins);
    xlabel('Synaptic weight bins');
    ylabel('Frequency');
    title('Synaptic Weight Distribution');

    subplot(2,3,[2 6]);
    raster(pt_output_spikes,testPeriodMS,DT);
    saveas(gcf,strcat('PT',simrun),'epsc');
    eps2pdf(strcat('PT',simrun,'.eps'),gspath);
    saveas(gcf,strcat('PT',simrun),'png');
end

% Plot training stimuli
for e=0:loops-1
    figure();
    trnSpikes = dlmread(strcat('E',int2str(e),'L0ExcitSpikes.dat'));
    trnSpikes(:,1) = []; % Delete the first column
    raster(trnSpikes,trainPeriodMS,DT);
    saveas(gcf,strcat('E',int2str(e),'inputs',simrun),'epsc');
    eps2pdf(strcat('E',int2str(e),'inputs',simrun,'.eps'),gspath);
    saveas(gcf,strcat('E',int2str(e),'inputs',simrun),'png');
end

% Plot post-training results

figure();
raster(ESpikes{1},testPeriodMS,DT); %L0ExcitSpikes
saveas(gcf,strcat('inputs',simrun),'epsc');
eps2pdf(strcat('inputs',simrun,'.eps'),gspath);
saveas(gcf,strcat('inputs',simrun),'png');

if plotInhib == 1
    if inputInhib == 1
        figure();
        raster(ISpikes{1},testPeriodMS,DT); %L0InhibSpikes
        saveas(gcf,strcat('InhibInputs',simrun),'epsc');
        eps2pdf(strcat('InhibInputs',simrun,'.eps'),gspath);
        saveas(gcf,strcat('InhibInputs',simrun),'png');
    end

    % Plot output layer inhibitory spikes
    figure();
    raster(ISpikes{nLayers},testPeriodMS,DT); %eval(['L',int2str(nLayers-1),'InhibSpikes'])
    saveas(gcf,strcat('InhibOutputs',simrun),'epsc');
    eps2pdf(strcat('InhibOutputs',simrun,'.eps'),gspath);
    saveas(gcf,strcat('InhibOutputs',simrun),'png');
end

figure();
subplot(2,3,1);
%imagesc(sort_matrix(L1weightsEfE'));
imagesc(eval(['L',int2str(nLayers-1),'weightsEfE']));
xlabel('Input neuron');
ylabel('Output neuron'); % Check these labels
title('Synaptic Weight Matrix');

subplot(2,3,4);
hist(reshape(eval(['L',int2str(nLayers-1),'weightsEfE']),1,[]),wbins);
xlabel('Synaptic weight bins');
ylabel('Frequency');
title('Synaptic Weight Distribution');

subplot(2,3,[2 6]);
raster(ESpikes{nLayers},testPeriodMS,DT); %eval(['L',int2str(nLayers-1),'ExcitSpikes'])
saveas(gcf, simrun, 'epsc');
eps2pdf(strcat(simrun,'.eps'),gspath);
saveas(gcf, simrun, 'png');



%% Calculate MD5 hashs
%eval('! md5 bweights.dat | tee bweights.md5');
%eval('! md5 output_spikes.dat | tee output_spikes.md5');
%eval('! md5 inhibitory_spikes.dat | tee inhibitory_spikes.md5');
%% Calculate mean and SD for each weight series for each presynaptic neuron for each loop and output record
% RWSTATS{NWLAYERS,NRECORDS_PL}=0;
% for lay=1:NWLAYERS
%     for r=1:NRECORDS_PL
%         wstats=zeros(NEXCIT,2,LOOPS);
%         for l=1:LOOPS
%             wstats(:,1,l)=mean(RWEIGHTS{lay,r}(:,:,l),2);
%             wstats(:,2,l)=std(RWEIGHTS{lay,r}(:,:,l),0,2);
%         end
%         RWSTATS{lay,r} = wstats;
%     end
% end
%% Clean up
clear remain
clear token

save matlab_workspace

%%
% for f=1:NRECORDS_PL
%     cell_analysis(f,1);
% end


function fRates = calc_fRates(spikeTimes, nStimuli, nTransPS, nCells, transTS)

fRates = zeros(nTransPS, nCells, nStimuli);
% spikeTimes = ESpikes{nLayers}; %eval(strcat('L',int2str(sLayer),'ExcitSpikes'));

edges = (0 : (nTransPS * nStimuli)) * transTS; % Calculate bin edges
spikeTimes(spikeTimes==0) = -1; % Don't count zeros padding the matrix
freqs = histc(spikeTimes',edges);
freqs(end-1,:)=freqs(end-1,:)+freqs(end,:); % include spikes at last timestep
freqs = freqs(1:end-1,:); % remove n^th +1 bin (spikeTimes==edges(end))
%freqs = freqs * 1/transP_Test; %-timewin makes this redundant

for s=1:nStimuli
    fRates(:,:,s) = freqs(1+(s-1)*nTransPS:s*nTransPS,:);
end
