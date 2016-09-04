function SpikeAnalysis(rdir,simrun)

tic;
if ~strcmp(rdir(end),filesep); rdir = strcat(rdir,filesep); end

if nargin == 1
    %rdir = '~/Results/SpikeNet/';
    content = dir(rdir);
    tmpcell = struct2cell(content);
    [~,ind] = sort_nat(tmpcell(1,:));
    clear tmpcell;
    for d = 1:length(content)
        simrun = content(ind(d)).name;
        if (content(ind(d)).isdir)
            ignoreDir = any(strcmpi(simrun, {'private','CVS','.','..'}));
            ignorePrefix = any(strncmp(simrun, {'@','.'}, 1));
            if ~(ignoreDir || ignorePrefix)
                SpikeAnalysis(rdir,simrun);
            end
        end
    end
    return;
end

if isa(simrun, 'numeric'); simrun = int2str(simrun); end

disp(['*** Evaluating: ',rdir,simrun,' ***']);
eval(['cd ',strcat(rdir,simrun)]);
gspath='/usr/local/bin/gs';

parameters;
if exist('test.m','file'); test; MP.printConnections=0; end % Legacy parameters files
if exist('DT','var'); MP.DT=DT; end
%if (useFilteredImages); system('rsync -L ~/bin/SpikeNet/wdir/image*.m .'); imageParams; end

close all;

%% Load all spike time data

%%% Analysis options %%%
plotInhib=false;
plotInput=true;
%plotOutput=true;
plotConnections=false;
plotTraining=true;
plotRecords=false;
savePDF=false; % Use exportfig() (in ~/Matlab/)
%%%%%%%%%%%%%%%%%%%%%%%%

% Compare print() to saveas() http://cvlab.epfl.ch/~ksmith/tips.html
% To make a Matlab axis expand to fill an entire figure window so the gray
% border is removed:
% >> set(gca, 'Position', [0 0 1 1]);
% For an image, you can call imshow with the following option:
% >> imshow(<image>, 'Border', 'tight');

%if (useFilteredImages); plotInput=false; end
outLayer = MP.nLayers - 1;
% Place into cells?
if MP.pretrain == 1
    if exist('preTraining.tbz','file')
        system('tar -xvf preTraining.tbz'); 
    end
    ptESpikes = cell(1,MP.nLayers);
    ptWeights = cell(1,MP.nLayers-1);
    ptFRates = cell(1,MP.nLayers);
    if (MP.SOM || strcmp(MP.axonDelay,'SOM'))
        ptSOMFRates = cell(1,MP.nLayers);
    end
end
if exist('postTraining.tbz','file')
    system('tar -xvf postTraining.tbz'); 
end
ESpikes = cell(1,MP.nLayers);
ISpikes = cell(1,MP.nLayers);
Weights = cell(1,MP.nLayers-1);
FRates = cell(1,MP.nLayers);
if (MP.SOM || strcmp(MP.axonDelay,'SOM'))
    SOMFRates = cell(1,MP.nLayers);
end
if MP.train==true
    if exist('training.tbz','file')
        system('tar -xvf training.tbz'); 
    end
    trnESpikes = cell(MP.loops,MP.nLayers);
    for e=0:MP.loops-1
        for l=0:MP.nLayers-1
            trnESpikes{e+1,l+1} = dlmread(strcat('E',int2str(e),'L',int2str(l),'ExcitSpikes.dat'));
            trnESpikes{e+1,l+1}(:,1) = []; % Delete the first column
        end
    end
end
if MP.printConnections && exist('connectivity.tbz','file')
    system('tar -xvf connectivity.tbz');
    AffEI = cell(1,MP.nLayers);
    AffElE = cell(1,MP.nLayers);
    AffIE = cell(1,MP.nLayers);
    AffII = cell(1,MP.nLayers);
    AffEfE = cell(1,MP.nLayers-1); % Make consistent with Weights cell?
end
for l=0:outLayer
    prefix = strcat(['L',int2str(l)]);
    if MP.printConnections && exist('connectivity.tbz','file')
        AffEI{l+1} = dlmread(strcat(prefix,'affNeuronsEI.dat')); %eval(['load ',prefix,'affNeuronsEI.dat']);
        AffElE{l+1} = dlmread(strcat(prefix,'affNeuronsElE.dat')); %eval(['load ',prefix,'affNeuronsElE.dat']);
        AffIE{l+1} = dlmread(strcat(prefix,'affNeuronsIE.dat')); %eval(['load ',prefix,'affNeuronsIE.dat']);
        AffII{l+1} = dlmread(strcat(prefix,'affNeuronsII.dat')); %eval(['load ',prefix,'affNeuronsII.dat']);
        if (l>0)
            AffEfE{l} = dlmread(strcat(prefix,'affNeuronsEfE.dat')); %eval(['load ',prefix,'affNeuronsEfE.dat']);
        end
    end
    if (l>0)
        Weights{l} = dlmread(strcat(prefix,'weightsEfE.dat')); %eval(['load ',prefix,'weightsEfE.dat']);
        Weights{l}(Weights{l}==0) = NaN;
        if MP.pretrain == 1
            ptWeights{l} = dlmread(strcat('pt',prefix,'weightsEfE.dat'));
            ptWeights{l}(ptWeights{l}==0) = NaN;
        end
    end
    if MP.pretrain == 1
        ptESpikes{l+1} = dlmread(strcat('pt',prefix,'ExcitSpikes.dat'));
        ptESpikes{l+1}(:,1) = []; % Delete first column
    end
    ESpikes{l+1} = dlmread(strcat(prefix,'ExcitSpikes.dat'));
    ESpikes{l+1}(:,1) = []; % Delete first column
    ISpikes{l+1} = dlmread(strcat(prefix,'InhibSpikes.dat'));
    ISpikes{l+1}(:,1) = []; % Delete first column
end

MP.layDim = cell(1,MP.nLayers);
for l=1:MP.nLayers
    if (MP.vSquare(l)==1)
        MP.layDim{l} = [round(sqrt(MP.vExcit(l))),round(sqrt(MP.vExcit(l)))];
    else
        MP.layDim{l} = [1,MP.vExcit(l)];
    end 
end

%% Plot connectivity
if MP.printConnections && plotConnections
%     cnxfig = figure();
%     maximize(cnxfig);
    


    ticks = 0:0.25:1;
    
%     % 3D plot
%     % Pass figure handle and plot all connections to/from a particular neuron

    
    % Plot EfE connections with preTraining weights
    cnxfig = figure();
    maximize(cnxfig);
    %subplot(3,1,1);
    nColours = plotConnectivity(AffEfE,ptWeights,MP); %[cnxfig,nColours]
    colorbar('YTick',1+ceil(nColours*ticks),'YTickLabels',ticks);
    title('Excitatory feed-forward pre-training weights');
    saveas(cnxfig, strcat('connectPreTrain',simrun),'png');
    close(cnxfig);
    
    % Plot EfE connections
    cnxfig = figure();
    maximize(cnxfig);
    %subplot(3,1,2);
    nColours = plotConnectivity(AffEfE,Weights,MP);
    %maximize(cnxfig);
    colorbar('YTick',1+ceil(nColours*ticks),'YTickLabels',ticks);
    title('Excitatory feed-forward post-training weights');
    saveas(cnxfig, strcat('connectPostTrain',simrun),'png');
    close(cnxfig);
    
    % Plot EfE weight changes
    cnxfig = figure();
    maximize(cnxfig);
    %subplot(3,1,3);
    deltaWeights = cell(1,MP.nLayers-1);
    for l=1:MP.nLayers-1
        deltaWeights{l} = (Weights{l} - ptWeights{l} + 1)/2; % Rescale: [0,1]
    end
    nColours = plotConnectivity(AffEfE,deltaWeights,MP);
    %maximize(cnxfig);
    colorbar('YTick',1+ceil(nColours*ticks),'YTickLabels',2*ticks-1);
    title('Excitatory feed-forward weight changes');
    saveas(cnxfig,strcat('connectWeightChange',simrun),'png');
    
end

%% Analyse firing rates
maxRate=200; % Last bin edge
step=10; % Spikes/s
thresh = 30; % Spikes/s for calculating transform overlaps

if MP.pretrain; ptSparseness = zeros(MP.nStimuli,MP.nTransPS,MP.nLayers); end
sparseness = zeros(MP.nStimuli,MP.nTransPS,MP.nLayers);
means = zeros(MP.nStimuli,MP.nTransPS,MP.nLayers);
stdds = zeros(MP.nStimuli,MP.nTransPS,MP.nLayers);
transTS = MP.transP_Test/MP.DT;

for l=1:MP.nLayers
    % Bin the Excitatory output spikes
    if MP.pretrain == 1
        ptCounts = calc_fRates(ptESpikes{l}, MP.nStimuli, MP.nTransPS, MP.vExcit(l), transTS);
        ptFRates{l} = ptCounts/MP.transP_Test; % Spikes/s
%        ptNRates = ptRates*refract;% Normalize according to maximum possible number of spikes
        
        if (MP.SOM || strcmp(MP.axonDelay,'SOM'))
            if (MP.vSquare(l))
                nCols = sqrt(MP.vExcit(l));
                nRows = nCols;
            else
                nCols = MP.vExcit(l);
                nRows = 1;
            end
            %SOMrates = permute(rates,[2,3,1]); % Necessary?
            %SOMrates = reshape(SOMrates, [nCols,nRows,MP.nTransPS,MP.nStimuli]);
            % Each layer of cells is transposed (C: row major, MATLAB: col major)

            ptSOMFRates{l} = zeros(nRows,nCols,MP.nTransPS,MP.nStimuli);
            for st=1:MP.nStimuli
                for tr=1:MP.nTransPS
                    ptSOMFRates{l}(:,:,tr,st) = reshape(ptFRates{l}(tr,:,st),[nCols,nRows])';
                end
            end


            fRates = figure();
            %mtit('PreTraining Firing Rates');
            maximize(fRates);
            %set(gca, 'Position', [0 0 1 1]);
            for st=1:MP.nStimuli
                for tr=1:MP.nTransPS
                    subplot(MP.nStimuli,MP.nTransPS,tr+(st-1)*MP.nTransPS);
                    %imagesc(SOMrates(:,:,tr,st)');
                    imagesc(ptSOMFRates{l}(:,:,tr,st)');
                    set(gca,'YDir','normal'); % Normal y-axis
                    if (MP.vSquare(l)); axis square; end
                    set(gca,'xtick',[],'ytick',[]);
                    title(['L',int2str(l-1),'S',int2str(st),'T',int2str(tr)]);
                end
            end
            saveas(fRates,strcat('ptL',int2str(l-1),'SOMfRates'),'png'); %gcf
            close(fRates);
        end
    end
    counts = calc_fRates(ESpikes{l}, MP.nStimuli, MP.nTransPS, MP.vExcit(l), transTS);
    FRates{l} = counts/MP.transP_Test; % Spikes/s
%    nRates = rates*refract;% Normalize according to maximum possible number of spikes

    if (MP.SOM || strcmp(MP.axonDelay,'SOM'))
        if (MP.vSquare(l))
            nCols = sqrt(MP.vExcit(l));
            nRows = nCols;
        else
            nCols = MP.vExcit(l);
            nRows = 1;
        end
        %SOMrates = permute(rates,[2,3,1]); % Necessary?
        %SOMrates = reshape(SOMrates, [nCols,nRows,MP.nTransPS,MP.nStimuli]);
        % Each layer of cells is transposed (C: row major, MATLAB: col major)
        
        SOMFRates{l} = zeros(nRows,nCols,MP.nTransPS,MP.nStimuli);
        for st=1:MP.nStimuli
            for tr=1:MP.nTransPS
                SOMFRates{l}(:,:,tr,st) = reshape(FRates{l}(tr,:,st),[nCols,nRows])';
            end
        end


        SOMfRates = figure();
        %mtit('PostTraining Firing Rates');
        maximize(SOMfRates);
        %set(gca, 'Position', [0 0 1 1]);
        for st=1:MP.nStimuli
            for tr=1:MP.nTransPS
                subplot(MP.nStimuli,MP.nTransPS,tr+(st-1)*MP.nTransPS);
                %imagesc(SOMrates(:,:,tr,st)');
                imagesc(SOMFRates{l}(:,:,tr,st)');
                set(gca,'YDir','normal'); % Normal y-axis
                if (MP.vSquare(l)); axis square; end
                set(gca,'xtick',[],'ytick',[]);
                title(['L',int2str(l-1),'S',int2str(st),'T',int2str(tr)]);
            end
        end
        saveas(SOMfRates,strcat('L',int2str(l-1),'SOMfRates'),'png');
        close(SOMfRates);
    end
    
    %set(gca,'position',get(0,'screensize'));
    %set(0,'units','normalized');
    %figure('Position',[0 0 1 1]);
    %set(0,'defaultfigureposition',[0 0 1 1]); 
    
    fRdist = figure();
    maximize(fRdist);
    %set(gca, 'Position', [0 0 1 1]);
    mtit(['L',int2str(l-1),' Firing Rate distributions']);
    disp(['L',int2str(l-1),': Overlap for successive transforms.']);
    for st=1:MP.nStimuli
        for tr=1:MP.nTransPS
            if MP.pretrain == 1
                % Calculate sparseness
                ptSTRates = ptFRates{l}(tr,:,st); % Single transform rates
                ptSparseness(st,tr,l) = (sum(ptSTRates)/MP.vExcit(l))^2 / (ptSTRates*ptSTRates' /MP.vExcit(l));
            end
            
            % Calculate overlap
            if tr>1
                last = find(FRates{l}(tr-1,:,st)>=thresh);
                next = find(FRates{l}(tr,:,st)>=thresh);
                overlap = intersect(last,next);
                disp(['S',int2str(st),'T',int2str(tr-1),' --> ','S',int2str(st),'T',int2str(tr),': ',...
                    num2str(100*length(overlap)/length(last)),'% overlap.']);
            end
            
            % Calculate sparseness
            stRates = FRates{l}(tr,:,st); % Single transform rates
            means(st,tr,l) = mean(stRates);
            stdds(st,tr,l) = std(stRates);
            sparseness(st,tr,l) = (sum(stRates)/MP.vExcit(l))^2 / (stRates*stRates' /MP.vExcit(l));
            
            edges = 0:step:maxRate+step;
            edges(end)=99999;%inf;
            n=histc(FRates{l}(tr,:,st),edges);
            edges(end)=maxRate+step;
            subplot(MP.nStimuli,MP.nTransPS,tr+(st-1)*MP.nTransPS);
            if ~MP.useFilteredImages
                bar(edges,n);
            else % Just plot all but the first bin and display #&% silent cells
                bar(edges(2:end),n(2:end)); 
            end
            %max2 = max(n(n~=max(n)));
            %breakplot(edges,n,35,round(0.99*MP.vExcit(l))-5,'Line',100.0);
            %breakplot(edges,n,1.1*max2,round(0.99*MP.vExcit(l)),'RPatch',80);
            title({['L',int2str(l-1),'S',int2str(st),'T',int2str(tr),' Silent:',num2str(100*n(1)/MP.vExcit(l),3),'%']; ...
                ['a:',num2str(sparseness(st,tr,l),3),' mean:',num2str(means(st,tr,l),3),' sd:',num2str(stdds(st,tr,l),3)]}); %' Tot:',int2str(sum(counts(tr,:,st)))
            xlim([-step/2,step/2 + maxRate]);
            %ylim([0,MP.vExcit(l)]);
            xlabel('Spikes/s (min edge)');
            %set(gca,'ytick',[]);
            %if tr==1
            %    set(gca, 'ytickMode','auto');
            %    ylabel('Frequency');
            %end


            %             [row,col] = find(ESpikes{l}>=start & ESpikes{l}<start+transTS & ESpikes{l}~=0);
            %             totSpikes = length(sub2ind(size(ESpikes{l}),row,col));
            %             %totSpikes = length(find(ESpikes{l}>=start & ESpikes{l}<start+transTS & ESpikes{l}~=0));
            %             counts=zeros(size(ESpikes{l},1));
            %             for n=1:vExcit(l)
            %                 counts(n)=length(find(ESpikes{l}(n,:)>=start & ESpikes{l}(n,:)<start+transTS & ESpikes{l}(n,:)~=0));
            %             end
        end
    end
    saveas(fRdist,strcat('L',int2str(l-1),'Distribution'),'png');
    close(fRdist);
end


if MP.pretrain
    disp('preTraining Sparseness layer-by-layer');
    disp(ptSparseness); 
end
disp('postTraining Sparseness layer-by-layer');
disp(sparseness);


%% Information theory analysis

% bpath = '~/bin/analysis/Release/';
sLayer = MP.nLayers-1; % Selected layer equal to output layer in 'C' numbering
% maxBins = nTransPS;
% minTrans = 4; %nTransPS;
% maxTrans = nTransPS;
% decoding = 1;
nBestCells = 5;

% Bin the Excitatory output spikes
if MP.pretrain == 1
    ptoutFrates = calc_fRates(ptESpikes{MP.nLayers}, MP.nStimuli, MP.nTransPS, MP.vExcit(MP.nLayers), MP.transP_Test/MP.DT); %* 1/DT
    ptINFO = infoAnalysis(ptoutFrates, MP.nStimuli, MP.nTransPS, MP.vExcit(MP.nLayers), nBestCells, sLayer, MP.transP_Test, simrun);
end
outFrates = calc_fRates(ESpikes{MP.nLayers}, MP.nStimuli, MP.nTransPS, MP.vExcit(MP.nLayers), MP.transP_Test/MP.DT);
INFO = infoAnalysis(outFrates, MP.nStimuli, MP.nTransPS, MP.vExcit(MP.nLayers), nBestCells, sLayer, MP.transP_Test, simrun);
 

%% Plot results of information analysis

plotMutual=false;
% Mutual single cell information
% plot(sort(msci{4},'descend')); %Sort 4th column (corrected info) and plot
if (plotMutual)
    mutInfo = figure();
    plot(sort(INFO{1},'descend')); %Sort 4th column (corrected info) and plot
    if MP.pretrain == 1
        hold on
        plot(sort(ptINFO{1},'descend'),'--r');
        legend('Trained','Untrained');
        hold off
    end
    axis([1 MP.vExcit(MP.nLayers) 0 1.05*log2(MP.nStimuli)]);
    xlabel('Cell rank');
    ylabel('Information (bits)');
    title('Mutual single cell information');
    fname = strcat('infoMutual',simrun);
    if (savePDF)
        saveas(mutInfo,fname,'epsc');
        eps2pdf(strcat(fname,'.eps'),gspath);
    end
    saveas(mutInfo,fname,'png');
end

% Stimulus specific single cell information
% plot(sort(sssci{2},'descend')); %Sort 2nd column and plot
singInfo = figure();
plot(sort(INFO{2},'descend')); %Sort 2nd column and plot
if MP.pretrain == 1
    hold on
    plot(sort(ptINFO{2},'descend'),'--r');
    legend('Trained','Untrained');
    hold off
end
axis([1 MP.vExcit(MP.nLayers) 0 1.05*log2(MP.nStimuli)]);
xlabel('Cell rank');
ylabel('Information (bits)');
title('Stimulus specific single cell information');
fname = strcat('infoSingle',simrun);
if (savePDF)
    saveas(singInfo,fname,'epsc');
    eps2pdf(strcat(infoSingle,'.eps'),gspath);
end
saveas(singInfo,fname,'png');

% Multiple cell information
% plot(mci{3}); %Plot 3rd column (already sorted ascending)
multInfo = figure();
plot(INFO{3}(2:end)); %Plot 3rd column (already sorted ascending)
if MP.pretrain == 1
    hold on
    plot(ptINFO{3}(2:end),'--r'); % First entry is 0 corresponding to 0 cells
    legend('Trained','Untrained');
    hold off
end
axis([1 (nBestCells*MP.nStimuli) 0 1.05*log2(MP.nStimuli)]);
xlabel('Number of best cells');
ylabel('Information (bits)');
title('Multiple cell information');
fname = strcat('infoMultiple',simrun);
if (savePDF)
    saveas(multInfo,fname,'epsc');
    eps2pdf(strcat(fname,'.eps'),gspath);
end
saveas(multInfo,fname,'png');

%% Perform network analysis

sampInt = 0.002; % bin width in seconds

testPeriodMS = MP.transP_Test * MP.nStimuli * MP.nTransPS * 1000;
trainPeriodMS = MP.transP_Train * MP.nStimuli * MP.nTransPS * 1000;
testPeriod = testPeriodMS/1000; % seconds
testBins = sampInt/2 : sampInt : (testPeriod - sampInt/2); % specify centres for hist(), or edges for histc()
Fs = 1/sampInt; % Sampling frequency (Hz)
NyquistLim = Fs/2; % (1/sampInt)/2
disp(['Nyquist Limit = ',num2str(NyquistLim,3),' Hz']);

for l=1:MP.nLayers
    spikes = ESpikes{l}.*MP.DT; % Spike times in seconds
    spikes(spikes(:)==0)=[];
    spikes=spikes(:)'; % Converts matrix to vector % reshape(A, 1, prod(size(A)))
    
    freqA = figure();
    
    subplot(2,1,1);
    freqs = hist(spikes,testBins);
    %freqs = histc(spikes,sampInt*(0:1:testPeriodMS/1000-1));
    %freqs = freqs - mean(freqs);
    NFFT = 2^nextpow2(length(freqs)); % Next power of 2 from length of freqs
    FFT = fft(freqs,NFFT);%/length(freqs);%Fs;
    %FFT(1) = []; % This is simply the sum of the data and can be removed
    % Plot first half of FFT (second half is a reflection of the first)
    f = NyquistLim * linspace(0,1,NFFT/2 + 1); % Gives f up to Nyquist limit NyquistLim * (0:NFFT-1); %
    plot(f,2*abs(FFT(1:NFFT/2 +1))/NFFT);       ylabel('|Amplitute|');
    %plot(f,(abs(FFT(1:NFFT/2 + 1)).^2)/NFFT);   ylabel('Power');
    xlabel('Frequency (Hz)');
    xlim([0,NyquistLim]); %xlim([-1,NyquistLim+1]);
    title(['L',int2str(l-1),' Fourier Spectrum']);
    
    subplot(2,1,2);    
    freqs = freqs / (MP.vExcit(l)*sampInt);
    %freqs = hist(spikes,testBins)/(MP.vExcit(l)*sampInt); % Normalise
    %freqs = freqs - mean(freqs);
    bar(testBins*1000,freqs); % Plot as msec
    xlim([0,testPeriodMS+(step/2)]);
    xlabel('Time (ms)'); %Seconds (Binned)
    %ylabel('Average Spikes/s'); % Without normalising would just be frequency
    ylabel(['L',int2str(l-1),'Total Spikes/s']);
    title(['L',int2str(l-1),' Excitatory Spike Histogram']);
    fname = strcat('L',int2str(l-1),'ExcitPopulation',simrun);
    if (savePDF)
        saveas(freqA,fname,'epsc');
        eps2pdf(strcat(fname,'.eps'),gspath);
    end
    saveas(freqA,fname,'png');
end

%% Plot Rasters and weights

% Calculate weight bins and centres
%wbins = 0.025:0.05:0.975;

% Plot pre-training results
if MP.pretrain==1
    preAnalysis = figure();
    subplot(2,3,1);
    imagesc(ptWeights{MP.nLayers-1}); 
    set(gca,'YDir','normal'); % Normal y-axis
    xlabel('Input neuron');
    ylabel('Output neuron'); % Check these labels
    title('Synaptic Weight Matrix');

    subplot(2,3,4);
    weights = reshape(ptWeights{MP.nLayers-1},1,[]);
    hist(weights,floor(sqrt(length(weights)))); % hist(reshape(ptWeights{MP.nLayers-1},1,[]),wbins);
    xlabel('Synaptic weight bins');
    ylabel('Frequency');
    title('Synaptic Weight Distribution');

    subplot(2,3,[2 6]);
    raster(ptESpikes{MP.nLayers},testPeriodMS,MP.DT); %pt_output_spikes
    fname = strcat('ptOutputAnalysis',simrun);
    if (savePDF)
        saveas(preAnalysis,fname,'epsc');
        eps2pdf(strcat(fname,'.eps'),gspath);
    end
    saveas(preAnalysis,fname,'png');
end
    
% Plot training stimuli (first layer only)
if plotTraining==true && MP.train==true
    for e=0:MP.loops-1
        trainInputs = figure();
        maximize(trainInputs);
        raster(trnESpikes{e+1,1},trainPeriodMS,MP.DT); % Only plotting inputs
        grid('minor');
        if (savePDF)
            saveas(trainInputs,strcat('E',int2str(e),'inputs',simrun),'epsc');
            eps2pdf(strcat('E',int2str(e),'inputs',simrun,'.eps'),gspath);
        end
        saveas(trainInputs,strcat('E',int2str(e),'inputs',simrun),'png');
        if e<MP.loops-1; close(trainInputs); end % Leave the last one
    end
end

% Plot post-training results
if plotInput==true
    postInputs = figure();
    raster(ESpikes{1},testPeriodMS,MP.DT); %L0ExcitSpikes
    if (savePDF)
        saveas(postInputs,strcat('inputs',simrun),'epsc');
        eps2pdf(strcat('inputs',simrun,'.eps'),gspath);
    end
    saveas(postInputs,strcat('inputs',simrun),'png');
    close(postInputs);
end

if plotInhib == 1
    if MP.inputInhib == 1
        inpInhib = figure();
        raster(ISpikes{1},testPeriodMS,MP.DT); %L0InhibSpikes
        if (savePDF)
            saveas(inpInhib,strcat('InhibInputs',simrun),'epsc');
            eps2pdf(strcat('InhibInputs',simrun,'.eps'),gspath);
        end
        saveas(inpInhib,strcat('InhibInputs',simrun),'png');
        close(inpInhib);
    end

    % Plot output layer inhibitory spikes
    outInhib = figure();
    raster(ISpikes{MP.nLayers},testPeriodMS,MP.DT); %eval(['L',int2str(nLayers-1),'InhibSpikes'])
    if (savePDF)
        saveas(outInhib,strcat('InhibOutputs',simrun),'epsc');
        eps2pdf(strcat('InhibOutputs',simrun,'.eps'),gspath);
    end
    saveas(outInhib,strcat('InhibOutputs',simrun),'png');
end

postAnalysis = figure();
subplot(2,3,1);
imagesc(Weights{MP.nLayers-1});
set(gca,'YDir','normal'); % Normal y-axis
xlabel('Input neuron');
ylabel('Output neuron'); % Check these labels
title('Synaptic Weight Matrix');

subplot(2,3,4);
weights = reshape(Weights{MP.nLayers-1},1,[]);
hist(weights,floor(sqrt(length(weights)))); %hist(reshape(Weights{MP.nLayers-1},1,[]),wbins);
xlabel('Synaptic weight bins');
ylabel('Frequency');
title('Synaptic Weight Distribution');

subplot(2,3,[2 6]);
raster(ESpikes{MP.nLayers},testPeriodMS,MP.DT); %eval(['L',int2str(nLayers-1),'ExcitSpikes'])
fname = strcat('OutputAnalysis',simrun);
if (savePDF)
    saveas(postAnalysis, fname, 'epsc');
    eps2pdf(strcat(fname,'.eps'),gspath);
end
saveas(postAnalysis, fname, 'png');

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
clear remain token

if (savePDF)
    system('mkdir eps');
    system('mv *.eps eps');
    system('mkdir pdf');
    system('mv *.pdf pdf');
end
system('mkdir infoAnalysis');
system('mv infos.* infom.* res* mr pcc sd_results beststiminfo.dat infoAnalysis');
% if (MP.nRecordsPL)
%     system('mkdir Records');
%     system('mv R*.dat Records');
% end
system('rm *.dat');



%% Cell Analysis

layer=outLayer; % C number
regime = 'Training';
loop=0;
%ticks = 0:0.25:1;
    
if (MP.nRecordsPL && plotRecords)
    records;
    system('mkdir Records');
    if exist('records.tbz','file')
        system('tar -xvf records.tbz'); 
    end
    system('mv R*.dat Records');
    %for l=0:MP.nLayers-1
    for r=1:MP.nRecordsPL
        neuron = MP.Records{c2m(layer)}(r);
        if (strcmpi(regime, 'preTesting'))
            spikes = ptESpikes{c2m(layer)}(c2m(neuron),:);
            prefix = strcat('Records/Rpt');
            weights = ptWeights; % Post training
        elseif (strcmpi(regime, 'Training'))
            spikes = trnESpikes{c2m(loop),c2m(layer)}(c2m(neuron),:);
            prefix = strcat(['Records/RE',int2str(loop)]);
            % Plot evolution of weights video.
            weights = Weights; % Tempory placeholder...
        elseif (strcmpi(regime, 'Testing'))
            spikes = ESpikes{c2m(layer)}(c2m(neuron),:);
            prefix = strcat('Records/R');
            weights = Weights; % Post training
        else
            disp('Error: Unknown regime!');
        end
        cellAnalysis(MP,layer,neuron,spikes,prefix,loop,MP.DT,TotalMS);
    
        % 3D plot - Pass figure handle and plot all connections to/from a particular neuron
        cnx3D = figure();
        plotNet3D(MP, weights, AffEfE, cnx3D, layer, neuron); % Check this is neuron ID not record ID
        % Change camera angle
        saveas(cnx3D, strcat('connectivity3D',simrun),'png'); % Post training
        %close(cnx3D);
    end
    %end
    plotNet3D(MP, weights, AffEfE); % Full network plot
    system('rm *.dat');
end

save matlab_workspace
disp(['Results saved in ', '<a href = "file://',pwd,filesep,'matlab_workspace.mat">',strcat(rdir,simrun),filesep,'matlab_workspace.mat</a>']);
toc;

% Matlab <--> C index conversions to make the code clearer and remind me!
% Extend to deal with row major (C) <--> column major (Matlab)?
function index = c2m(index)

index=index+1;

function index = m2c(index)

index=index-1;

function fRates = calc_fRates(spikeTimes, nStimuli, nTransPS, nCells, transTS)

fRates = zeros(nTransPS, nCells, nStimuli);
edges = (0 : (nTransPS * nStimuli)) * transTS; % Calculate bin edges
if isvector(spikeTimes) % Pad into matrix
    spikeTimes = [spikeTimes, zeros(length(spikeTimes),1)]; % This forces histc to bin spikes for each cell even when there is at most one spike per cell
end
spikeTimes(spikeTimes==0) = -1; % Don't count zeros padding the matrix
freqs = histc(spikeTimes',edges); %n(k) counts the value x(i) if edges(k) <= x(i) < edges(k+1). 
%The last bin counts any values of x that match edges(end).
%freqs(end-1,:)=freqs(end-1,:)+freqs(end,:); % include spikes at last
%timestep - not necessary since using C numbering last possible spike is
%edges(end)-1
freqs(end,:) = []; % remove n^th +1 bin (spikeTimes==edges(end)) %freqs = freqs(1:end-1,:); 
%freqs = freqs * 1/transP_Test; %-timewin makes this redundant
for s=1:nStimuli
    fRates(:,:,s) = freqs(1+(s-1)*nTransPS:s*nTransPS,:);
end

%return fRates;


function nColours = plotConnectivity(cnx, weights, MP)

LWS = 4; % Line width scale
mSize = 5;

%handle = figure();
colormap(jet);
cmp = colormap;
nColours=size(cmp,1);
scaleWidth=false;

for l=0:MP.nLayers-1
    % Plot excitatory neurons
    hold on;
    plot(0:MP.vExcit(l+1)-1,l,'ok','LineStyle','none','MarkerFaceColor','g','MarkerSize',mSize); % use draw?
    
    if (l>0)
        for n=1:MP.vExcit(l+1)
            nSyn = cnx{l}(n,1);
            for s=1:nSyn
                line = plot([n-1,cnx{l}(n,s+1)],[l,l-1],'Color',cmp(ceil(weights{l}(n,s)*(nColours)),:),'LineWidth',LWS/2); %,'LineWidth',LWS*weights{l}(n,s)); % ([Xdata],[Ydata])
                if (scaleWidth)
                    set(line,'LineWidth', LWS * weights{l}(n,s));
                end
            end
        end
    end
    %colorbar;
    set(gca,'YTick',0:MP.nLayers-1);
    xlabel('Neuron');
    ylabel('Layer');
    box on;
    hold off;
end

axis([-0.5 max(MP.vExcit)-0.5 -0.5 MP.nLayers-0.5]);
%return handle;


function [xcoords,ycoords] = calcNetCoords(layDim) %, markProp)
%net=figure();
%nodes=[4,5]; markProp='og';
sx=1/layDim(2); % nCols
sy=1/layDim(1); % nRows
xcoords=(sx/2)+(sx*(0:m2c(layDim(2))));
%ycoords=(sy/2)+(0:layDim(1)-1)*sy;
ycoords=(sy/2)+(sy*(m2c(layDim(1)):-1:0)); % Order is reversed
if length(ycoords)==1 % Row vector layer
    ycoords=ycoords*ones(length(xcoords));
elseif length(xcoords)==1
    xcoords=xcoords*ones(length(ycoords));
end
%zcoords=layer*ones(length(xcoords));

%[X,Y,Z] = meshgrid(xcoords, ycoords, zcoords);

%plot(repmat(xcoords,nodes(2),1),repmat(ycoords,nodes(1),1)','og')
%plot(X,Y,markProp);
%xlim([0,1]);
%ylim([0,1]);
%     sx=1/layDim{1}(2); % nCols
%     sy=1/layDim{1}(1); % nRows
%     x1=(sx/2)+(0:layDim{1}(2)-1)*sx;
%     y1=(sy/2)+(0:layDim{1}(1)-1)*sy;
%     y1=y1*ones(length(x1));
%     z1=0*ones(1,length(x1));

% narray[l][n].x = (sp_x/2 + (narray[l][n].col * sp_x)) * mp->spatialScale; // col: 0,1,...,rowLen-1
% narray[l][n].y = (sp_y/2 + (((colLen - 1) - narray[l][n].row) * sp_y)) * mp->spatialScale; // row indices and y coordinates run opposite


function fh = plotNet3D(MP, weights, cnx, fh, layer, neurons)
% Use C indexing i.e. pass layer and neurons starting at 0
% 3D plot: Plot all connections to/from a particular neuron
if nargin <= 3
    % Plot whole net
    fh = figure();
    wholeNet = true;
    recurse = false;
else
%     if length(neuron)>1
%         % [row,col]
%     end
    figure(fh);
    wholeNet = false;
    recurse = true;
end

ticks = 0:0.25:1;
hold on
colormap(jet);
cmp = colormap;
nColours=size(cmp,1);

%weights = Weights; % Post training


coords = cell(MP.nLayers,2);%3);

% Calculate coordinates and plot neurons  % Redo indexing!!!!!
for l=1:MP.nLayers
    [x,y] = calcNetCoords(MP.layDim{l});
    coords{l,1} = x * MP.spatialScale;
    coords{l,2} = y * MP.spatialScale;
    [X,Y,Z] = meshgrid(coords{l,1},coords{l,2},m2c(l));
    plot3(X,Y,Z,'ok');
    if wholeNet
        if l>1
            neurons = 0:MP.vExcit(l)-1;
            plotPreCnx3D(MP, fh, cmp, nColours, coords, cnx, weights, m2c(l), neurons, recurse);
%             for n=1:MP.vExcit(l) % post synaptic
%                 [c,r] = ind2sub(MP.layDim{l},n); % swap indices to transpose (c is row major)
%                 nSyn = cnx{l-1}(n,1);
% %                 disp(['n=',int2str(n),', r=',int2str(r),', c=',int2str(c),', nSyn=',int2str(nSyn)]);
%                 for s=1:nSyn
%                     preN = c2m(cnx{l-1}(n,s+1));
%                     [preC,preR] = ind2sub(MP.layDim{l-1},preN); % swap indices to transpose (c is row major)
%                     line([coords{l,1}(c),coords{l-1,1}(preC)],[coords{l,2}(r),coords{l-1,2}(preR)],[l-1,l-2],'Color',cmp(ceil(weights{l-1}(n,s)*(nColours)),:)); % Change s
%                 end
%             end
        end
    end
end

if ~wholeNet
    if layer>0  % Draw all presynaptic connections
        plotPreCnx3D(MP, fh, cmp, nColours, coords, cnx, weights, layer, neurons, recurse);
    end
    if layer<MP.nLayers-1  % Draw all postsynaptic connections
        
    end
end

hold off
colorbar('YTick',1+ceil(nColours*ticks),'YTickLabels',ticks);
xlabel('x'); xlim([0,MP.spatialScale]);
ylabel('y'); ylim([0,MP.spatialScale]);
zlabel('Layer'); zlim([0,m2c(MP.nLayers)]);

% Outside
% title('Excitatory feed-forward post-training weights');
% saveas(cnx3D, 'connectivity3D','png'); % Post training
% close(cnx3D);

function plotPreCnx3D(MP, fh, cmp, nColours, coords, cnx, weights, l, neurons, recurse)

figure(fh);
for n=1:length(neurons)
    [c,r] = ind2sub(MP.layDim{c2m(l)},c2m(neurons(n)));
    nSyn = cnx{l}(c2m(neurons(n)),1);
    preNeurons = cnx{l}(c2m(neurons(n)),1+(1:nSyn)); % List of connected presyn neurons
    for s=1:nSyn %length(preNeurons)
        preN = c2m(preNeurons(s));
        [preC,preR] = ind2sub(MP.layDim{l},preN); % swap indices to transpose (c is row major)
        line([coords{c2m(l),1}(c),coords{c2m(l-1),1}(preC)],[coords{c2m(l),2}(r),coords{c2m(l-1),2}(preR)],[l,l-1],'Color',cmp(ceil(weights{l}(c2m(neurons(n)),s)*(nColours)),:)); % Change s
    end
    if l>1 && recurse==true
        plotPreCnx3D(MP, fh, cmp, nColours, coords, cnx, weights, l-1, preNeurons); % cnx{layer-1}(1+(1:nSyn))
    end
        
end

%%
%function animateWeights()
% 3 options
% 1) Save frames for later playback as a movie
% e.g.
% for k = 1:16
%  plot(fft(eye(k+16)))
%  axis equal
%  M(k) = getframe;
% end
% movie(M,30) % Playback 30 times

% 2) Continually erase/redraw 

% 3) Redefine XData, YData, ZData, and/or CData plot object properties use refreshdata
% First call set to define a data source for an axis, e.g.,
% set(obj_handle,'YDataSource','varname') to make a persistent link
% Then call refreshdata after you or your code updates varname.

% Change camera angle (and rotate throughout movie)

function orbit(deg)
[az el] = view;
rotvec = 0:deg/10:deg;
for i = 1:length(rotvec)
    view([az+rotvec(i) el])
    drawnow
end
