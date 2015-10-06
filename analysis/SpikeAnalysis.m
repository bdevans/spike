function fh = SpikeAnalysis(rdir,simrun,reprocess,plotLevel)
if nargin == 0 %# If data is not specified return a cell array of function handles
    fh = {@c2m @m2c @binSpikes @plotFourierSpectrum @psth @calc_fRates @plotConnectivity @calcNetCoords @plotNet3D @plotPreCnx3D};
    return %# Return from the function
end

if nargin < 4
    plotLevel = 1; % Save PNGs and standard plots only
    if nargin < 3
        reprocess=false;
    end
end

global useCorrFunc; % Check this: http://www.mathworks.co.uk/help/techdoc/ref/global.html
[success,~]=license('checkout','Econometrics_Toolbox');
if (success) % autocorr, crosscorr
    useCorrFunc = true;
else % Use signal_toolbox (xcorr)
    useCorrFunc = false;
end

global useStatsFunc; %normcdf in info analysis routines
[success,~]=license('checkout','Statistics_Toolbox'); 
if (success); useStatsFunc = true; else useStatsFunc = false; end

tic;
if ~strcmp(rdir(end),filesep); rdir = strcat(rdir,filesep); end
%prefix = '~/model/results/';
%if length(rdir)<=length(prefix) || ~strcmpi(prefix,rdir(1:length(prefix))); rdir=['~/model/results/',rdir]; end

if nargin == 1 || strcmpi(simrun,'all') % Pass directory path
    ignoreDirs = {'private','CVS','.','..','plots'}; %Directories to ignore
    content = dir(rdir);
    tmpcell = struct2cell(content);
    [~,ind] = sort_nat(tmpcell(1,:));
    clear tmpcell;
    for d = 1:length(content)
        simrun = content(ind(d)).name;
        if any(strncmp(simrun, {'@','.'}, 1)) % Files to ignore (redundant)
            continue;
        end
        if (content(ind(d)).isdir && ~any(strcmpi(simrun, ignoreDirs)))
            SpikeAnalysis(rdir,simrun,reprocess,plotLevel);
        end
    end
    if plotLevel
        genReport(rdir);
    end
    return;
end

if (isa(simrun, 'numeric') && length(simrun) > 1) % Pass a vector of sims
    for e=1:length(simrun)
        SpikeAnalysis(rdir,simrun(e),reprocess,plotLevel);
    end
    if plotLevel
        genReport(rdir);
    end
    return;
end

%% Main routine begins

if isa(simrun, 'numeric'); simrun = int2str(simrun); end


if exist([rdir,simrun,filesep,'matlab_workspace.mat'],'file') && ~reprocess
    disp(['*** Skipping: ',rdir,simrun,' - already processed! ***']);
    return;
else
    disp(['*** Evaluating: ',rdir,simrun,' ***']);
end

cd(strcat(rdir,simrun)); 

INF = Inf;  % Boilerplate for infinite condSpeed

% Load Model parameters
parameters;

% Legacy parameter fixes
if (exist('DT','var') && ~isfield(MP,'DT')); MP.DT=DT; end
%if (useFilteredImages); system('rsync -L ~/bin/SpikeNet/wdir/image*.m .'); imageParams; end
if ~isfield(MP,'nTestStimuli'); MP.nTestStimuli = MP.nStimuli; end
if ~isfield(MP,'nTestTransPS'); MP.nTestTransPS = MP.nTransPS; end
if ~isfield(MP,'stimGroups'); MP.stimGroups = false; end
if ~isfield(MP,'priorPhases'); MP.priorPhases = false; end
if ~isfield(MP,'trainElE'); MP.trainElE = false; end
if ~exist('TestMS','var'); TestMS = MP.transP_Test * MP.nTestStimuli * MP.nTestTransPS * 1000; end
if ~exist('EpochMS','var'); EpochMS = MP.transP_Train * MP.nStimuli * MP.nTransPS * 1000; end
if (MP.initElE==3); MP.initElE = 'SOM'; end
if (MP.initEfE==3); MP.initEfE = 'SOM'; end
if any(strcmpi({MP.initEfE,MP.initElE,MP.axonDelay},'SOM')); MP.SOM = true; end
if strcmpi(MP.axonDelay,'MinD'); MP.axonDelay = 0; end

if plotLevel > 0
    close all;
    if exist('CLIparams.m','file')
        popupmessage('CLIparams.m',['CLI Parameters: ',rdir,' [',simrun,']']);
    else
        popupmessage('parameters.m',['Parameters: ',rdir,' [',simrun,']']);
    end
end

%% Analysis options %%
% 0: No plots (for regenerating mat files)
% 1: Standard plots & save pngs (& report)
% 2: Extended plots & save pdfs/pngs
% 3: Full plots & save pdfs/pngs

savePDF=false;
plotInhib=false;
plotInput=false;
plotConnections=false;
plotTraining=false;
plotFullTraining=false;
plotRecords=false;
plotArchitecture=false;
plotMutual=false;

if plotLevel == 0
    disp('Skipping plots!');
    
elseif plotLevel >= 1
    plotInput=true;
    plotTraining=true;
    plotRecords=true;
    
    if plotLevel >= 2
       savePDF=true;
       plotInhib=true;
       plotConnections=true;

        if plotLevel >= 3
            plotFullTraining=true;
            plotArchitecture=true;
            plotMutual=true;
            if plotLevel > 3
                disp(['Undefined plotLevel: ',int2str(plotLevel)]);
                return;
            end
        end
    end
else
    disp(['Undefined plotLevel: ',int2str(plotLevel)]);
    return;
end
    

if (MP.useFilteredImages); plotTraining=false; end %plotInput
%%%%%%%%%%%%%%%%%%%%%%%%

% Replace epsc wih epsc2 in saveas commands (also strcat with [])
% Try laprint for latex: http://www.uni-kassel.de/fb16/rat/matlab/laprint/
% exportfig() is used at fMRIB: http://www.mathworks.co.uk/company/newsletters/digest/june00/export/index.html http://www.mathworks.co.uk/company/newsletters/digest/december00/export.html
%  e.g. exportfig(gcf,'test.eps','color','rgb','resolution',300)
% Try export_fig and savefig
% http://pgfplots.sourceforge.net/
% http://www.mathworks.com/matlabcentral/fileexchange/21286-matlabfrag
% http://win.ua.ac.be/~nschloe/content/matlab2tikz
% http://andrewjpage.com/index.php?/archives/27-Using-an-EPS-in-Latex.html
% Tips: http://blogs.mathworks.com/loren/2007/12/11/making-pretty-graphs/

% Compare print() to saveas() http://cvlab.epfl.ch/~ksmith/tips.html

% For an image, you can call imshow with the following option:
% >> imshow(<image>, 'Border', 'tight');

% To expand a figure to fill an entire window (remove gray border & labels):
% >> set(gca, 'Position', [0 0 1 1]);
% To crop the figure tight up to the labels and title (eps2pdf does this): 
% set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
% http://nibot-lab.livejournal.com/73290.html

    %set(gca,'position',get(0,'screensize'));
    %set(0,'units','normalized');
    %figure('Position',[0 0 1 1]);
    %set(0,'defaultfigureposition',[0 0 1 1]); 

%% Load all spike trains, weights and connectivity data

disp(['[',simrun,'] ','Loading Data']); % Add timestamps

outLayer = m2c(MP.nLayers); %MP.nLayers - 1;

if MP.pretrain == 1
    if exist('preTraining.tbz','file')
        system('tar -xvf preTraining.tbz'); 
    end
    ptESpikes = cell(1,MP.nLayers);
    ptWeightsEfE = cell(1,MP.nLayers-1);
    if MP.trainElE
        ptWeightsElE = cell(1,MP.nLayers);
    end
    ptFRates = cell(1,MP.nLayers);
    if (MP.SOM) % || MP.initElE==3 || strcmp(MP.axonDelay,'SOM'))
        ptSOMFRates = cell(1,MP.nLayers);
    end
end
if exist('postTraining.tbz','file')
    system('tar -xvf postTraining.tbz'); 
end
ESpikes = cell(1,MP.nLayers);
ISpikes = cell(1,MP.nLayers);
WeightsEfE = cell(1,MP.nLayers-1);
if MP.trainElE
    WeightsElE = cell(1,MP.nLayers);
end
FRates = cell(1,MP.nLayers);
if (MP.SOM) % || MP.initElE==3 || strcmp(MP.axonDelay,'SOM'))
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
            if sum(trnESpikes{e+1,l+1}(:,1)>0) % Legacy fix
            for n=1:MP.vExcit(c2m(l))
                trnESpikes{e+1,l+1}(n,trnESpikes{e+1,l+1}(n,1)+2:end) = -1; %NaN;
            end
            else
                trnESpikes{e+1,l+1}(trnESpikes{e+1,l+1}(:,:)==0) = -1; %NaN;
            end
            trnESpikes{e+1,l+1}(:,1) = []; % Delete the first column
        end
    end
end

if MP.priorPhases % Create cell arrays for ElE_ training and testing
    if exist('PP_stimuli.m','file')
        pStr = 'PP_';
        newDATnames = true;
    else
        pStr = '_PP_'; % Legacy
        newDATnames = false;
    end
    %if (oldPP); pStr = ['_',pStr]; end % Legacy

    if exist('PPpreTraining.tbz','file')
        system('tar -xvf PPpreTraining.tbz');
    end
    PP_WeightsElE = cell(1,MP.nLayers);
    if MP.pretrain
        ptPP_ESpikes = cell(1,MP.nLayers);
        ptPP_WeightsEfE = cell(1,MP.nLayers-1);
        ptPP_WeightsElE = cell(1,MP.nLayers);
        ptPP_FRates = cell(1,MP.nLayers);
        if (MP.SOM) % || MP.initElE==3 || strcmp(MP.axonDelay,'SOM'))
            ptPP_SOMFRates = cell(1,MP.nLayers);
        end
    end
    % Post training
    if exist('PPtraining.tbz','file')
        system('tar -xvf PPtraining.tbz');
    end
    PP_ESpikes = cell(1,MP.nLayers);
    PP_ISpikes = cell(1,MP.nLayers);
    PP_WeightsEfE = cell(1,MP.nLayers-1);
    PP_WeightsElE = cell(1,MP.nLayers);
    PP_FRates = cell(1,MP.nLayers);
    if (MP.SOM) % || MP.initElE==3 || strcmp(MP.axonDelay,'SOM'))
        PP_SOMFRates = cell(1,MP.nLayers);
    end
    if MP.train
%         if exist('training.tbz','file')
%             system('tar -xvf training.tbz');
%         end
        PP_trnESpikes = cell(MP.loops,MP.nLayers);
        for e=0:MP.loops-1
            for l=0:MP.nLayers-1
                PP_trnESpikes{e+1,l+1} = dlmread(strcat(pStr,'E',int2str(e),'L',int2str(l),'ExcitSpikes.dat'));
                if sum(PP_trnESpikes{e+1,l+1}(:,1)>0) % Legacy fix
                    for n=1:MP.vExcit(c2m(l))
                        PP_trnESpikes{e+1,l+1}(n,PP_trnESpikes{e+1,l+1}(n,1)+2:end) = -1; %NaN;
                    end
                else
                    PP_trnESpikes{e+1,l+1}(PP_trnESpikes{e+1,l+1}(:,:)==0) = -1; %NaN;
                end
                PP_trnESpikes{e+1,l+1}(:,1) = []; % Delete the first column
            end
        end
    end
end

if MP.printConnections % Load Connectivity
    if exist('connectivity.tbz','file')
        system('tar -xvf connectivity.tbz');
    end
    AffEI = cell(1,MP.nLayers);
    AffElE = cell(1,MP.nLayers);
    AffIE = cell(1,MP.nLayers);
    AffII = cell(1,MP.nLayers);
    AffEfE = cell(1,MP.nLayers-1); % Make consistent with Weights cell?
    
    if MP.axonDelay
        if MP.delayEfE
            AffDelaysEfE = cell(1,MP.nWLayers);
        end
        if MP.delayElE
            AffDelaysElE = cell(1,MP.nLayers);
        end
        if MP.delayEI
            AffDelaysEI = cell(1,MP.nLayers);
        end
    end
end

for l=0:outLayer
    prefix = strcat(['L',int2str(l)]);
    if MP.printConnections %&& exist('connectivity.tbz','file')
        AffEI{l+1} = dlmread(strcat(prefix,'affNeuronsEI.dat')); %eval(['load ',prefix,'affNeuronsEI.dat']);
        AffElE{l+1} = dlmread(strcat(prefix,'affNeuronsElE.dat')); %eval(['load ',prefix,'affNeuronsElE.dat']);
        AffIE{l+1} = dlmread(strcat(prefix,'affNeuronsIE.dat')); %eval(['load ',prefix,'affNeuronsIE.dat']);
        AffII{l+1} = dlmread(strcat(prefix,'affNeuronsII.dat')); %eval(['load ',prefix,'affNeuronsII.dat']);
        if (l>0)
            AffEfE{l} = dlmread(strcat(prefix,'affNeuronsEfE.dat')); %eval(['load ',prefix,'affNeuronsEfE.dat']);
            for n=1:MP.vExcit(c2m(l)) % Legacy fix required?
                AffEfE{l}(n,AffEfE{l}(n,1)+2:end) = -1;
            end
        end
        for n=1:MP.vExcit(c2m(l)) % Legacy fix required?
            AffElE{l+1}(n,AffElE{l+1}(n,1)+2:end) = -1;
            AffIE{l+1}(n,AffIE{l+1}(n,1)+2:end) = -1;
        end
        for n=1:MP.vInhib(c2m(l))
            AffEI{l+1}(n,AffEI{l+1}(n,1)+2:end) = -1; %NaN;
            AffII{l+1}(n,AffII{l+1}(n,1)+2:end) = -1;
        end
        % Do not clip the first column (nSyn) as 3D plots use it
        
        % Load distances? Upper triangle is zero padded
        %distElE{l+1} = dlmread(['L',int2str(l),'distElE.dat']); 
        
        if MP.axonDelay % So far delays are only on efferents of Excit neurons
            if ~strcmp(MP.delayEfE,'MinD') && l>0
                AffDelaysEfE{l} = dlmread([prefix,'affDelaysEfE.dat']);
                for n=1:MP.vExcit(c2m(l)) % Legacy fix required?
                    AffDelaysEfE{l}(n,AffDelaysEfE{l}(n,1)+2:end) = -1;
                end
            end
            if ~strcmp(MP.delayElE,'MinD')
                AffDelaysElE{l+1} = dlmread([prefix,'affDelaysElE.dat']);
                for n=1:MP.vExcit(c2m(l)) % Legacy fix required?
                    AffDelaysElE{l+1}(n,AffDelaysElE{l+1}(n,1)+2:end) = -1;
                end
            end
            if ~strcmp(MP.delayEI,'MinD')
                AffDelaysEI{l+1} = dlmread([prefix,'affDelaysEI.dat']);
                for n=1:MP.vInhib(c2m(l)) % Legacy fix required?
                    AffDelaysEI{l+1}(n,AffDelaysEI{l+1}(n,1)+2:end) = -1;
                end
            end
        end
    end
    
    if MP.priorPhases
        if exist('PPpostTraining.tbz','file')
            system('tar -xvf PPpostTraining.tbz');
        end
        if MP.pCnxElE(c2m(l)) > eps
            PP_WeightsElE{c2m(l)} = dlmread([pStr,prefix,'weightsElE.dat']);
            if sum(PP_WeightsElE{c2m(l)}(:,1)>0) % Legacy fix
                for n=1:MP.vExcit(c2m(l))
                    PP_WeightsElE{c2m(l)}(n,PP_WeightsElE{c2m(l)}(n,1)+2:end) = NaN;
                end
                PP_WeightsElE{c2m(l)}(:,1) = []; % Delete the first column
            else
                PP_WeightsElE{c2m(l)}(PP_WeightsElE{c2m(l)}==0) = NaN; %(:)==0
            end
            if MP.pretrain
                if newDATnames
                    ptPP_WeightsElE{c2m(l)} = dlmread([pStr,'pt',prefix,'weightsElE.dat']);
                else
                    ptPP_WeightsElE{c2m(l)} = dlmread(['pt',pStr,prefix,'weightsElE.dat']);
                end
                if sum(ptPP_WeightsElE{c2m(l)}(:,1)>0) % Legacy fix
                    for n=1:MP.vExcit(c2m(l))
                        ptPP_WeightsElE{c2m(l)}(n,ptPP_WeightsElE{c2m(l)}(n,1)+2:end) = NaN;
                    end
                    ptPP_WeightsElE{c2m(l)}(:,1) = []; % Delete the first column
                else ptPP_WeightsElE{c2m(l)}(ptPP_WeightsElE{c2m(l)}==0) = NaN;
                end
            end
        end
        if (l>0) % EfE weights
            PP_WeightsEfE{l} = dlmread(strcat(pStr,prefix,'weightsEfE.dat')); %eval(['load ',prefix,'weightsEfE.dat']);
            if sum(PP_WeightsEfE{l}(:,1)>0) % Legacy fix
                for n=1:MP.vExcit(c2m(l))
                    PP_WeightsEfE{l}(n,PP_WeightsEfE{l}(n,1)+2:end) = NaN;
                end
                PP_WeightsEfE{l}(:,1) = []; % Delete the first column
            else PP_WeightsEfE{l}(PP_WeightsEfE{l}==0) = NaN;
            end
            if MP.pretrain == 1
                if newDATnames
                    ptPP_WeightsEfE{l} = dlmread(strcat(pStr,'pt',prefix,'weightsEfE.dat'));
                else
                    ptPP_WeightsEfE{l} = dlmread(strcat('pt',pStr,prefix,'weightsEfE.dat'));
                end
                if sum(ptPP_WeightsEfE{l}(:,1)>0) % Legacy fix
                    for n=1:MP.vExcit(c2m(l))
                        ptPP_WeightsEfE{l}(n,ptPP_WeightsEfE{l}(n,1)+2:end) = NaN;
                    end
                    ptPP_WeightsEfE{l}(:,1) = []; % Delete the first column
                else ptPP_WeightsEfE{l}(ptPP_WeightsEfE{l}==0) = NaN;
                end
            end
        end
        
        if MP.pretrain % Load Spikes
            if newDATnames
                ptPP_ESpikes{l+1} = dlmread(strcat(pStr,'pt',prefix,'ExcitSpikes.dat'));
            else
                ptPP_ESpikes{l+1} = dlmread(strcat('pt',pStr,prefix,'ExcitSpikes.dat'));
            end
            if sum(ptPP_ESpikes{l+1}(:,1)>0) % Legacy fix
                for n=1:MP.vExcit(c2m(l))
                    ptPP_ESpikes{l+1}(n,ptPP_ESpikes{l+1}(n,1)+2:end) = -1;
                end
            else
                ptPP_ESpikes{l+1}(ptPP_ESpikes{l+1}(:,:)==0) = -1;
            end
            ptPP_ESpikes{l+1}(:,1) = []; % Delete first column
        end
        PP_ESpikes{l+1} = dlmread(strcat(pStr,prefix,'ExcitSpikes.dat'));
        if sum(PP_ESpikes{l+1}(:,1)>0) % Legacy fix
            for n=1:MP.vExcit(c2m(l))
                PP_ESpikes{l+1}(n,PP_ESpikes{l+1}(n,1)+2:end) = -1;
            end
        else
            PP_ESpikes{l+1}(PP_ESpikes{l+1}(:,:)==0) = -1;
        end
        PP_ESpikes{l+1}(:,1) = []; % Delete first column
        
        PP_ISpikes{l+1} = dlmread(strcat(pStr,prefix,'InhibSpikes.dat'));
        if sum(PP_ISpikes{l+1}(:,1)>0) % Legacy fix
            for n=1:MP.vInhib(c2m(l))
                PP_ISpikes{l+1}(n,PP_ISpikes{l+1}(n,1)+2:end) = -1;
            end
        else
            PP_ISpikes{l+1}(PP_ISpikes{l+1}(:,:)==0) = -1;
        end
        PP_ISpikes{l+1}(:,1) = []; % Delete first column
    end % End of loading ElE training and testing data
    
    
    if (MP.pCnxElE(c2m(l)) > eps) && MP.trainElE && ~MP.priorPhases % Reconsider.........................
        WeightsElE{c2m(l)} = dlmread([prefix,'weightsElE.dat']);
        if sum(WeightsElE{c2m(l)}(:,1)>0) % Legacy fix
            for n=1:MP.vExcit(c2m(l))
                WeightsElE{c2m(l)}(n,WeightsElE{c2m(l)}(n,1)+2:end) = NaN;
            end
            WeightsElE{c2m(l)}(:,1) = []; % Delete the first column
        else WeightsElE{c2m(l)}(WeightsElE{c2m(l)}==0) = NaN; %(:)==0
        end
        if MP.pretrain
            ptWeightsElE{c2m(l)} = dlmread(['pt',prefix,'weightsElE.dat']);
            if sum(ptWeightsElE{c2m(l)}(:,1)>0) % Legacy fix
                for n=1:MP.vExcit(c2m(l))
                    ptWeightsElE{c2m(l)}(n,ptWeightsElE{c2m(l)}(n,1)+2:end) = NaN;
                end
                ptWeightsElE{c2m(l)}(:,1) = []; % Delete the first column
            else ptWeightsElE{c2m(l)}(ptWeightsElE{c2m(l)}==0) = NaN;
            end
        end
    end
    
    if (l>0) % EfE weights
        WeightsEfE{l} = dlmread(strcat(prefix,'weightsEfE.dat')); %eval(['load ',prefix,'weightsEfE.dat']);
        if sum(WeightsEfE{l}(:,1)>0) % Legacy fix
            for n=1:MP.vExcit(c2m(l))
                WeightsEfE{l}(n,WeightsEfE{l}(n,1)+2:end) = NaN;
            end
            WeightsEfE{l}(:,1) = []; % Delete the first column
        else WeightsEfE{l}(WeightsEfE{l}==0) = NaN;
        end
        if MP.pretrain == 1
            ptWeightsEfE{l} = dlmread(strcat('pt',prefix,'weightsEfE.dat'));
            if sum(ptWeightsEfE{l}(:,1)>0) % Legacy fix
                for n=1:MP.vExcit(c2m(l))
                    ptWeightsEfE{l}(n,ptWeightsEfE{l}(n,1)+2:end) = NaN;
                end
                ptWeightsEfE{l}(:,1) = []; % Delete the first column
            else ptWeightsEfE{l}(ptWeightsEfE{l}==0) = NaN;
            end
        end
    end
    
    if MP.pretrain == 1
        ptESpikes{l+1} = dlmread(strcat('pt',prefix,'ExcitSpikes.dat'));
        if sum(ptESpikes{l+1}(:,1)>0) % Legacy fix
        for n=1:MP.vExcit(c2m(l))
            ptESpikes{l+1}(n,ptESpikes{l+1}(n,1)+2:end) = -1;
        end
        else ptESpikes{l+1}(ptESpikes{l+1}(:,:)==0) = -1;
        end
        ptESpikes{l+1}(:,1) = []; % Delete first column
    end
    
    ESpikes{l+1} = dlmread(strcat(prefix,'ExcitSpikes.dat'));
    if sum(ESpikes{l+1}(:,1)>0) % Legacy fix
    for n=1:MP.vExcit(c2m(l))
        ESpikes{l+1}(n,ESpikes{l+1}(n,1)+2:end) = -1;
    end
    else ESpikes{l+1}(ESpikes{l+1}(:,:)==0) = -1;
    end
    ESpikes{l+1}(:,1) = []; % Delete first column
    
    if l>0 || MP.inputInhib
        ISpikes{l+1} = dlmread(strcat(prefix,'InhibSpikes.dat'));
        if sum(ISpikes{l+1}(:,1)>0) % Legacy fix
            for n=1:MP.vInhib(c2m(l))
                ISpikes{l+1}(n,ISpikes{l+1}(n,1)+2:end) = -1;
            end
        else ISpikes{l+1}(ISpikes{l+1}(:,:)==0) = -1;
        end
        ISpikes{l+1}(:,1) = []; % Delete first column
    end
end

if exist('stimuli.tbz','file')
    system('tar -xvf stimuli.tbz');
end
if exist('stimuli.m','file'); stimuli; end
%if (MP.K == MP.nTestStimuli); MP.nStimuli = 1; end % Legacy Hack

if (exist('STIM','var') && ~MP.useFilteredImages && isfield(STIM,'train'))
    incInj = true;
    if isfield(MP,'sInputs'); sInputs = MP.sInputs; else sInputs = MP.vExcit(1); end
    STIM.protoTest = zeros(MP.nTestStimuli,sInputs);
%     if MP.newTestSet
%         pats = STIM.test;
%     else
%         STIM.test = STIM.train;
%         pats = STIM.test; %STIM.train;
%     end
    if ~(MP.newTestSet || isfield('STIM','test'))
        STIM.test = STIM.train;
    end
    pats = STIM.test;
    for s=1:MP.nTestStimuli
        for t=1:MP.nTestTransPS
            STIM.protoTest(s,:) = or(STIM.protoTest(s,:), pats{s,t});
        end
    end
    
    STIM.protoTrain = zeros(MP.nStimuli,sInputs);
    for s=1:MP.nStimuli
        for t=1:MP.nTransPS
            STIM.protoTrain(s,:) = or(STIM.protoTrain(s,:), STIM.train{s,t});
        end
    end
else
    incInj = false;
end


if MP.stimGroups    % Load stimulus group prototypes
    if exist('prototypes.dat','file'); system('mv prototypes.{dat,stm}'); end % Legacy fix
    prototypes = load('prototypes.stm');
    gMembers = cell(MP.nGroups,1);
    for g=1:MP.nGroups
        gMembers{g} = find(prototypes(g,:));
    end
    gShared = find(sum(prototypes)>1);
end

if MP.K > 1 % Simultaneously presented training stimuli %change to sMembers? See below (stimSpikes)...
    sMembers = cell(MP.nTestStimuli,1);
    if MP.localRep % This assumes stimuli are a contiguous block of 1's
        block = floor(MP.nFiringNeurons + ((MP.nTransPS - 1) * MP.shift));
        for g=1:MP.nTestStimuli
            sMembers{g} = (g-1)*block+(1:block);
        end
%         for p=0:MP.nTestStimuli-1
%             for trans=0:MP.nTransPS-1
%                 for n=(p*block)+(trans*MP.shift):(MP.nFiringNeurons+(p*block)+(trans*MP.shift))-1
%                     tst_stimuli(p,trans,n) = MP.current;
%                 end
%             end
%         end
    end
end

MP.layDim = cell(1,MP.nLayers);
for l=1:MP.nLayers
    if (MP.vSquare(l)==1)
        MP.layDim{l} = [round(sqrt(MP.vExcit(l))),round(sqrt(MP.vExcit(l)))];
    else
        MP.layDim{l} = [1,MP.vExcit(l)];
    end 
end

if (MP.nRecords && ~exist('records.tbz','file'))% && ~exist('Records','dir'))
    mkdir('Records'); % Load records once and delete dat files
    [status,message,messageID] = movefile('R*.dat','Records/');
    if status == 0
        disp(message); disp(messageID);
    end
end
    
%system('rm *.dat');
% The following is safer, particularly for reanalysing old results
if exist('preTraining.tbz','file')
    system('tar -tf preTraining.tbz | xargs rm');
end
if exist('postTraining.tbz','file')
    system('tar -tf postTraining.tbz | xargs rm');
end
if exist('training.tbz','file')
    system('tar -tf training.tbz | xargs rm');
end
if exist('PPpreTraining.tbz','file')
    system('tar -tf PPpreTraining.tbz | xargs rm');
end
if exist('PPtraining.tbz','file')
    system('tar -tf PPtraining.tbz | xargs rm');
end
if exist('connectivity.tbz','file')
    system('tar -tf connectivity.tbz | xargs rm');
end
if exist('PPpostTraining.tbz','file')
    system('tar -tf PPpostTraining.tbz | xargs rm');
end
if exist('stimuli.tbz','file')
    system('tar -tf stimuli.tbz | xargs rm');
end

disp(['[',simrun,'] ','Data loaded']);


%% Plot 2D connectivity
if plotLevel && MP.printConnections && plotConnections  
    % Change to plot only recorded neurons with spatial scaling...
    ticks = 0:0.25:1;
    
    % Plot EfE connections with preTraining weights
    disp(['[',simrun,'] ','Plotting connectivity with pretraining EfE weights']);
    cnxfig = figure();    maximize(cnxfig);
    nColours = plotConnectivity(AffEfE,ptWeightsEfE,MP); %[cnxfig,nColours]
    colorbar('YTick',1+ceil(nColours*ticks),'YTickLabel',ticks);
    title('Excitatory feed-forward pre-training weights');
    saveFig(cnxfig,['connectPreTrain',simrun],plotLevel);    close(cnxfig);
    
    % Plot EfE connections
    disp(['[',simrun,'] ','Plotting connectivity with post-training EfE weights']);
    cnxfig = figure();    maximize(cnxfig);
    nColours = plotConnectivity(AffEfE,WeightsEfE,MP);
    colorbar('YTick',1+ceil(nColours*ticks),'YTickLabel',ticks);
    title('Excitatory feed-forward post-training weights');
    saveFig(cnxfig,['connectPostTrain',simrun],plotLevel);   close(cnxfig);
    
    % Plot EfE weight changes
    disp(['[',simrun,'] ','Plotting connectivity with EfE weight changes']);
    cnxfig = figure();    maximize(cnxfig);
    deltaWeights = cell(1,MP.nLayers-1);
    for l=1:MP.nLayers-1
        deltaWeights{l} = (WeightsEfE{l} - ptWeightsEfE{l} + 1)/2; % Rescale: [0,1]
    end
    nColours = plotConnectivity(AffEfE,deltaWeights,MP);
    colorbar('YTick',1+ceil(nColours*ticks),'YTickLabel',2*ticks-1);
    title('Excitatory feed-forward weight changes');
    saveFig(cnxfig,['connectWeightDiff',simrun],plotLevel);  close(cnxfig);
    
    disp(['[',simrun,'] ','Connectivity plotted successfully']);
end


%% Analyse firing rates
thresh  =    30; % Spikes/s for calculating transform overlaps

% Change to floor(sqrt(MP.vExcit(l))) for number of bins?
if MP.pretrain; ptSparseness = zeros(MP.nTestStimuli,MP.nTestTransPS,MP.nLayers); end
sparseness = zeros(MP.nTestStimuli,MP.nTestTransPS,MP.nLayers);
means = zeros(MP.nTestStimuli,MP.nTestTransPS,MP.nLayers);
stdds = zeros(MP.nTestStimuli,MP.nTestTransPS,MP.nLayers);
transTS = MP.transP_Test/MP.DT;
% s_N:'standard deviation of the sample' (uncorrected)
% s:  'sample standard deviation' (Bessel's correction - denominator (N-1))
% Default (or flag 0 not flag 1): s (corrected)

for l=1:MP.nLayers
    disp(['[',simrun,'] ','Layer ',int2str(l),' firing rates']);
    % Bin the Excitatory output spikes
    if MP.pretrain == 1
        ptCounts = calc_fRates(ptESpikes{l}, MP.nTestStimuli, MP.nTestTransPS, MP.vExcit(l), transTS);
        ptFRates{l} = ptCounts/MP.transP_Test; % Spikes/s
%        ptNRates = ptRates*refract;% Normalize according to maximum possible number of spikes
        
        if (MP.SOM) % || MP.initElE==3 || strcmp(MP.axonDelay,'SOM'))
            disp(['[',simrun,'] ','Analysing Pretraining SOM firing rates']);
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

            ptSOMFRates{l} = zeros(nRows,nCols,MP.nTestTransPS,MP.nTestStimuli);
            for st=1:MP.nTestStimuli
                for tr=1:MP.nTestTransPS
                    ptSOMFRates{l}(:,:,tr,st) = reshape(ptFRates{l}(tr,:,st),[nCols,nRows])';
                end
            end

            if plotLevel > 0
                fRates = figure('Name', 'Pre-training SOM firing rates');
                maximize(fRates); %set(gca, 'Position', [0 0 1 1]);
                maxFR = max(max(max(max(ptSOMFRates{l}(:,:,:,:)))));
                if (maxFR > 0)
                    for st=1:MP.nTestStimuli
                        for tr=1:MP.nTestTransPS
                            sh=subplot(MP.nTestStimuli,MP.nTestTransPS,tr+(st-1)*MP.nTestTransPS);
                            imagesc(ptSOMFRates{l}(:,:,tr,st)',[0 maxFR]); % Normalise colours
                            set(gca,'YDir','normal'); % Normal y-axis
                            set(gca,'xtick',[],'ytick',[]);
                            if (MP.vSquare(l));
                                axis square;
                            else
                                if (tr==1); 
                                    if (l==1) && ~MP.useFilteredImages
                                        set(gca,'YTick',0:MP.nFiringNeurons:MP.vExcit(l));
                                    else
                                        set(gca,'YTickMode','auto');
                                    end
                                end
                            end
                            if tr==MP.nTestTransPS;
                                ax=get(sh,'position');
                                colorbar; %ch=colorbar;    ylabel(ch,'Spikes/s');
                                set(sh,'position',ax);
                            end
                            title(['S',int2str(st),'T',int2str(tr)]); %title({['L',int2str(l-1),'S',int2str(st)];['T',int2str(tr)]});
                        end
                    end
                    saveFig(fRates,['ptL',int2str(l-1),'SOMfRates'],plotLevel);
                else
                    disp(['[',simrun,'] ','Silent neurons: skipping pretraining SOM firing rate plots']);
                end
                close(fRates);
            end
            disp(['[',simrun,'] ','Pretraining done']);
        end
    end
    
    counts = calc_fRates(ESpikes{l}, MP.nTestStimuli, MP.nTestTransPS, MP.vExcit(l), transTS);
    FRates{l} = counts/MP.transP_Test; % Spikes/s
%    nRates = rates*refract;% Normalize according to maximum possible number of spikes

    if (MP.SOM) % || MP.initElE==3 || strcmp(MP.axonDelay,'SOM'))
        disp(['[',simrun,'] ','Analysing Post-training SOM firing rates']);
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
        
        SOMFRates{l} = zeros(nRows,nCols,MP.nTestTransPS,MP.nTestStimuli);
        for st=1:MP.nTestStimuli
            for tr=1:MP.nTestTransPS
                SOMFRates{l}(:,:,tr,st) = reshape(FRates{l}(tr,:,st),[nCols,nRows])';
            end
        end

        if plotLevel > 0
            SOMfRates = figure('Name', 'Post-training SOM firing rates');
            maximize(SOMfRates); %set(gca, 'Position', [0 0 1 1]);
            maxFR = max(max(max(max(SOMFRates{l}(:,:,:,:)))));
            if (maxFR > 0)
                for st=1:MP.nTestStimuli
                    for tr=1:MP.nTestTransPS
                        sh=subplot(MP.nTestStimuli,MP.nTestTransPS,tr+(st-1)*MP.nTestTransPS);
                        imagesc(SOMFRates{l}(:,:,tr,st)',[0 maxFR]);
                        set(gca,'YDir','normal'); % Normal y-axis
                        set(gca,'xtick',[],'ytick',[]);
                        if (MP.vSquare(l));
                            axis square;
                        else
                            if (tr==1); 
                                if (l==1) && ~MP.useFilteredImages
                                    set(gca,'YTick',0:MP.nFiringNeurons:MP.vExcit(l)); 
                                else
                                    set(gca,'YTickMode','auto');
                                end 
                            end
                        end
                        if tr==MP.nTestTransPS;
                            ax=get(sh,'position');
                            colorbar; %ch=colorbar;    ylabel(ch,'Spikes/s');
                            set(sh,'position',ax);
                        end
                        title(['S',int2str(st),'T',int2str(tr)]); %title({['L',int2str(l-1),'S',int2str(st)];['T',int2str(tr)]});
                    end
                end
                saveFig(SOMfRates,['L',int2str(l-1),'SOMfRates'],plotLevel);
            else
                disp(['[',simrun,'] ','Silent neurons: skipping SOM firing rate plots']);
            end
            close(SOMfRates);
        end
        disp(['[',simrun,'] ','Post-training done']);
    end
    
    disp(['[',simrun,'] ','Analysing firing rate distributions']);
    for st=1:MP.nTestStimuli
        for tr=1:MP.nTestTransPS
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
        end
    end
    
    if plotLevel > 0
        disp(['[',simrun,'] ','L',int2str(m2c(l)),': Overlap for successive transforms.']);
        if ~MP.useFilteredImages || (MP.useFilteredImages && l~=0)
            minfrpl = 0; %min(min(min(FRates{l}))); % Minimum firing rate per layer
        else
            minfrpl = 1; % Omit quiescent cells
        end
        maxfrpl = max(max(max(FRates{l}))); % Maximum firing rate per layer
        if maxfrpl > 0
            fRdist = figure('Name', ['L',int2str(m2c(l)),' Firing Rate distributions']);
            maximize(fRdist); %set(gca, 'Position', [0 0 1 1]);
            edges = linspace(minfrpl, maxfrpl, floor(sqrt(MP.vExcit(l))));
            for st=1:MP.nTestStimuli
                for tr=1:MP.nTestTransPS
                    subplot(MP.nTestStimuli,MP.nTestTransPS,tr+(st-1)*MP.nTestTransPS);
                    hist(FRates{l}(tr,:,st),edges);
                    
                    if ~MP.useFilteredImages || (MP.useFilteredImages && l~=0)
                        ylim([0,MP.vExcit(l)]);
                        if tr~=1; set(gca,'YTickLabel',[]); end
                        title(['S',int2str(st),'T',int2str(tr)]);
                    else % Exclude the first bin and print #&% silent cells
                        title({['S',int2str(st),'T',int2str(tr),' Silent:',num2str(100*n(1)/MP.vExcit(l),3),'\%']; ...
                            ['a:',num2str(sparseness(st,tr,l),3),'$ \bar{x}:$',num2str(means(st,tr,l),3),' s:',num2str(stdds(st,tr,l),3)]}...
                            ,'interpreter','latex');
                    end
                    xlim([-1, 10*ceil(maxfrpl/10)+1]); % +/- 1 for bar width
                    if st==MP.nTestStimuli
                        xlabel('Spikes/s (min edge)');
                    else
                        set(gca,'XTickLabel',[]);
                    end
                    
                end
            end
            saveFig(fRdist,['L',int2str(l-1),'Distribution',simrun],plotLevel);
            close(fRdist);
            disp(['[',simrun,'] ','Layer ',int2str(l),' firing rates distributions plotted']);
        else
            disp(['[',simrun,'] ','Layer ',int2str(l),': No cells fired!']);
        end
    end
end


if MP.pretrain
    disp(['[',simrun,'] ','preTraining Sparseness layer-by-layer']);
    disp(ptSparseness); 
end
disp(['[',simrun,'] ','postTraining Sparseness layer-by-layer']);
disp(sparseness);


%% Information theory analysis
if ~MP.stimGroups || MP.priorPhases
% bpath = '~/bin/analysis/Release/';
sLayer = MP.nLayers-1; % Selected layer equal to output layer in 'C' numbering
% maxBins = nTransPS;
% minTrans = 4; %nTransPS;
% maxTrans = nTransPS;
% decoding = 1;
nBestCells = 5; % Per stimulus

% if MP.stimGroups
%     objects = MP.nGroups;
%     transforms = MP.nStimuli;
% else
%     objects = MP.nStimuli;
%     transforms = MP.nTransPS;
% end

disp(['[',simrun,'] ','Calculating information analyses']);
% Bin the Excitatory output spikes
if MP.pretrain == 1
    ptoutFrates = ptFRates{c2m(sLayer)}*MP.transP_Test; %calc_fRates(ptESpikes{MP.nLayers}, MP.nStimuli, MP.nTransPS, MP.vExcit(MP.nLayers), MP.transP_Test/MP.DT); %* 1/DT
    ptINFO = infoAnalysis(ptoutFrates, MP.nTestStimuli, MP.nTestTransPS, MP.vExcit(MP.nLayers), nBestCells, sLayer, MP.transP_Test, simrun);
    system('mkdir ptInfoAnalysis');
    system('mv infos.* infom.* res* mr pcc sd_results beststiminfo.dat ptInfoAnalysis');
    [ptInfo,ptBC] = calcInfo(ptESpikes{c2m(sLayer)},1000*MP.transP_Test,MP.vExcit(c2m(sLayer)),MP.nTestStimuli,MP.nTestTransPS,MP.transP_Test/MP.DT,false);
end
outFrates = FRates{c2m(sLayer)}*MP.transP_Test; %calc_fRates(ESpikes{MP.nLayers}, MP.nStimuli, MP.nTransPS, MP.vExcit(MP.nLayers), MP.transP_Test/MP.DT);
[INFO, bestCells] = infoAnalysis(outFrates, MP.nTestStimuli, MP.nTestTransPS, MP.vExcit(MP.nLayers), nBestCells, sLayer, MP.transP_Test, simrun);

system('mkdir infoAnalysis');
system('mv infos.* infom.* res* mr pcc sd_results beststiminfo.dat infoAnalysis');

disp(['[',simrun,'] ','Processing new Information analysis']);
[Info,BC] = calcInfo(ESpikes{c2m(sLayer)},1000*MP.transP_Test,MP.vExcit(c2m(sLayer)),MP.nTestStimuli,MP.nTestTransPS,MP.transP_Test/MP.DT,plotLevel);
% if ~isequal(INFO,nInfo)
%     disp(INFO);
%     disp(nInfo);
% end
if ~isequal(bestCells, BC)
    disp(bestCells);
    disp(BC);
else
    if ~isequal(INFO{2},Info.Single)
        disp([INFO{2},Info.Single]);
    end
end
disp('Compare res.stimtemp with this:');
disp(Info.Single);

disp(['[',simrun,'] ','Information analyses calculated']);

%% Plot results of information analysis

if plotLevel > 0 && (MP.nTransPS >= 4)
    
    disp(['[',simrun,'] ','Plotting information analyses']);
    % Mutual single cell information
    if (plotMutual)
        mutInfo = figure('Name', 'Mutual single cell information');
        plot(sort(INFO{1},'descend')); %Sort msci{4} 4th column (corrected info) and plot
        if MP.pretrain == 1
            hold on
            plot(sort(ptINFO{1},'descend'),'--r');
            legend('Trained','Untrained');
            hold off
        end
        axis([1 MP.vExcit(MP.nLayers) 0 1.05*log2(MP.nTestStimuli)]);
        xlabel('Cell rank');        ylabel('Information (bits)');
        title('Mutual single cell information');
        saveFig(mutInfo,['infoMutual',simrun],plotLevel);
    end
    
    % Stimulus specific single cell information
    singInfo = figure('Name', 'Stimulus specific single cell information');
    plot(sort(INFO{2},'descend')); %Sort sssci{2} 2nd column and plot
    if MP.pretrain == 1
        hold on
        plot(sort(ptINFO{2},'descend'),'--r');
        legend('Trained','Untrained');
        hold off
    end
    axis([1 MP.vExcit(MP.nLayers) 0 1.05*log2(MP.nTestStimuli)]);
    %set(gca,'XTick',0:sqrt(MP.vExcit(MP.nLayers)):MP.vExcit(MP.nLayers));
    xlabel('Cell rank');    ylabel('Information (bits)');
    title('Stimulus specific single cell information');
    saveFig(singInfo,['infoSingle',simrun],plotLevel);
    
    % Multiple cell information
    multInfo = figure('Name', 'Multiple cell information');
    plot(INFO{3}(2:end)); %Plot mci{3} 3rd column (already sorted ascending)
    if MP.pretrain == 1
        hold on
        plot(ptINFO{3}(2:end),'--r'); % First entry is 0 corresponding to 0 cells
        legend('Trained','Untrained');
        hold off
    end
    axis([1 (nBestCells*MP.nTestStimuli) 0 1.05*log2(MP.nTestStimuli)]);
    xlabel('Number of best cells');    ylabel('Information (bits)');
    title('Multiple cell information');
    saveFig(multInfo,['infoMultiple',simrun],plotLevel);
    disp(['[',simrun,'] ','Information analyses plotted']);

    
    % Plot firing rates for the best (most informative) cells
    disp(['[',simrun,'] ','Plotting firing rates of the most informative output cells']);
    fRates = figure('Name', 'Best Cell Firing rates');
    %maximize(fRates);
    for st=1:MP.nTestStimuli
        for n=1:nBestCells
            subplot(MP.nTestStimuli, nBestCells, n+(st-1)*nBestCells);
            plot(FRates{c2m(outLayer)}(:,bestCells(n,st),st),'x-');
            hold on
            if MP.pretrain
                plot(ptFRates{c2m(outLayer)}(:,bestCells(n,st),st),'--');
                %legend('Trained','Untrained');
            end
            for st2=1:MP.nTestStimuli
                if (st ~= st2)
                    plot(FRates{c2m(outLayer)}(:,bestCells(n,st),st2),'+:r');
                end
            end
            hold off
            xlabel('Transform');
            xlim([1 MP.nTestTransPS]);
            ylabel('Firing Rate (Spikes/s)');
            title(['L',int2str(outLayer),'N',int2str(m2c(bestCells(n,st))),' S',int2str(st)]);
        end
    end
    saveFig(fRates,['L',int2str(l-1),'BestCellfRates',simrun],plotLevel);
    %close(fRates);
    disp(['[',simrun,'] ','Most informative output cells plotted']);
    
end

end


%% Perform network analysis

binsize = 0.005; %sampInt = 0.002; % bin width in seconds

% Filtering??? e.g. filtfilt or sshist(); sskernel();

if plotLevel > 0
    disp(['[',simrun,'] ','Plotting population analyses (PSTH and Fourier spectrum)']);
    tLength = MP.transP_Test * MP.nTestStimuli * MP.nTestTransPS;
    for l=1:MP.nLayers
        disp(['[',simrun,'] ','Layer ',int2str(l),' post-training population analysis']);
        spikes = ESpikes{l}.*MP.DT; % Spike times in seconds
        [freqs,bins] = binSpikes(spikes, binsize, tLength);
        
        freqA = figure('Name', ['L',int2str(m2c(l)),' post-training population analysis']);
        subplot(2,1,1); % Plot Fourier Spectrum
        plotFourierSpectrum(freqs,binsize,freqA);
        title(['L',int2str(l-1),' Fourier Spectrum']);
        
        subplot(2,1,2); % Plot PSTH
        psth(freqs,binsize,bins,tLength,MP.vExcit(l),freqA);
        title(['L',int2str(l-1),' Excitatory Spike Histogram']);
        
        saveFig(freqA,['L',int2str(l-1),'ExcitPopulation',simrun],plotLevel);
        disp(['[',simrun,'] ','Layer ',int2str(l),' post-training done']);
    end
end

if MP.stimGroups % Plot inputs during testing
    tLength = TestMS / 1000;
    
    if MP.pretrain
        disp(['[',simrun,'] Inputs pre-training groups analysis']);
        
        if MP.priorPhases
            % Be careful with PPtLength: using nStimuli not nTestStimuli
            % because PP test stimulus is the second phase training stimulus
            PPtLength = MP.nStimuli * MP.nTransPS * MP.transP_Test;
            ptPP_GroupSpikes = cell(MP.nGroups,1);
            for g=1:MP.nGroups % Bin Spikes
                spikes = ptPP_ESpikes{1}(gMembers{g},:).*MP.DT;
                [ptPP_GroupSpikes{g},bins] = binSpikes(spikes,binsize,PPtLength);
            end
        end
        
        ptGroupSpikes = cell(MP.nGroups,1);
        for g=1:MP.nGroups % Bin Spikes
            spikes = ptESpikes{1}(gMembers{g},:).*MP.DT;
            [ptGroupSpikes{g},bins] = binSpikes(spikes,binsize,tLength);
        end
        
        if plotLevel > 0
            ptTestG = figure('Name', 'Inputs pre-training groups analysis');
            for g=1:MP.nGroups
                % Plot Fourier Spectrum
                subplot(MP.nGroups,3,3*g-2);%2,1+(g-1)*2);
                plotFourierSpectrum(ptGroupSpikes{g},binsize,ptTestG)
                title(['G',int2str(g),' Inputs Fourier Spectrum']);
                % Plot PSTH
                subplot(MP.nGroups,3,3*g-1);%2,2+(g-1)*2);
                psth(ptGroupSpikes{g},binsize,bins,tLength,MP.vExcit(l),ptTestG);
                title(['G',int2str(g),' Inputs Excitatory Spike Histogram']);
                % Plot autocorrelations
                subplot(MP.nGroups,3,3*g)
                nLags = MP.transP_Test / binsize;
                if useCorrFunc
                    autocorr(ptGroupSpikes{g});%,nLags); % Should this use binned spikes?
                    set(gca, 'XTickMode', 'manual'); % Stops saveas recomputing XTicks
                    set(gca,'XTickLabel',get(gca,'XTick')*binsize*1000); %Show ms not lags
                else
                    [c_ww,lags]=xcorr(ptGroupSpikes{g},nLags,'coeff');
                    stem(lags(nLags+1:end)*binsize*1000,c_ww(nLags+1:end)); % +ve lags are enough
                    ylabel('Correlation coefficient');
                end
                xlabel('Lag (ms)');
                title(['G',int2str(g),' Inputs Autocorrelation']);
            end
            saveFig(ptTestG,['ptGroupInputAnalysis',simrun],plotLevel);
            close(ptTestG);
        end
    end
    
    disp(['[',simrun,'] Inputs post-training groups analysis']);
    
    if MP.priorPhases
        PPtLength = MP.nStimuli * MP.nTransPS * MP.transP_Test;
        PP_GroupSpikes = cell(MP.nGroups,1);
        for g=1:MP.nGroups % Bin Spikes
            spikes = PP_ESpikes{1}(gMembers{g},:).*MP.DT; %prototypes(g,:)>0
            [PP_GroupSpikes{g},bins] = binSpikes(spikes,binsize,PPtLength);
        end
    end
    
    GroupSpikes = cell(MP.nGroups,1);
    for g=1:MP.nGroups % Bin Spikes
        spikes = ESpikes{1}(gMembers{g},:).*MP.DT; %prototypes(g,:)>0
        [GroupSpikes{g},bins] = binSpikes(spikes,binsize,tLength);
    end
    
    if plotLevel > 0
        testG = figure('Name', 'Inputs post-training groups analysis');
        maximize(testG);
        for g=1:MP.nGroups
            % Plot Fourier Spectrum
            subplot(MP.nGroups,3,3*g-2);%2,1+(g-1)*2);
            plotFourierSpectrum(GroupSpikes{g},binsize,testG)
            title(['G',int2str(g),' Inputs Fourier Spectrum']);
            % Plot PSTH
            subplot(MP.nGroups,3,3*g-1);%2,2+(g-1)*2);
            psth(GroupSpikes{g},binsize,bins,tLength,MP.vExcit(l),testG);
            title(['G',int2str(g),' Inputs Excitatory Spike Histogram']);
            % Plot autocorrelations
            subplot(MP.nGroups,3,3*g)
            nLags = MP.transP_Test / binsize;
            if useCorrFunc
                autocorr(GroupSpikes{g});%,nLags); % Should this use binned spikes?
                set(gca, 'XTickMode', 'manual'); % Stops saveas recomputing XTicks
                set(gca,'XTickLabel',get(gca,'XTick')*binsize*1000); %Show ms not lags
            else
                [c_ww,lags]=xcorr(GroupSpikes{g},nLags,'coeff');
                stem(lags(nLags+1:end)*binsize*1000,c_ww(nLags+1:end)); % +ve lags are enough
                ylabel('Correlation coefficient');
            end
            xlabel('Lag (ms)');
            title(['G',int2str(g),' Inputs Autocorrelation']);
        end
        saveFig(testG,['GroupInputAnalysis',simrun],plotLevel);
        close(testG);
    end
end

% if MP.K > 1 % 
%     tLength = TrainMS / 1000;
%     epoch = 1; % Use first epoch only
%     
%     if MP.pretrain
%         disp(['[',simrun,'] Inputs combined training stimuli analysis']);
%         
%         GroupSpikes = cell(MP.nTestStimuli,1);
%         for g=1:MP.nTestStimuli
%             spikes = ESpikes{1}(gMembers{g},:).*MP.DT; %members = find(prototypes(g,:));
%             [GroupSpikes{g},bins] = binSpikes(spikes,binsize,tLength);
%         end
%     end
%     % See below...
% end
        
if MP.train==true && MP.K > 1 %plotTraining==true && 
    % Split input spike trains up according to stimulus (amalgamate transforms)
    tLength = EpochMS / 1000;
    l=1; % Set layer to input layer
    stimSpikes = cell(MP.loops,MP.nStimuli);
    combined = cell(1,MP.loops);
    if (MP.K == MP.nTestStimuli) % Hack
        nStim = MP.nTestStimuli; 
    else
        nStim = MP.nStimuli; 
    end % Hack
    if (plotFullTraining)
        lastEpoch=MP.loops;
    else
        lastEpoch=1;
    end
    for e=1:MP.loops%lastEpoch
        disp(['[',simrun,'] ','Input layer training population analysis and autocorrelations, epoch ',int2str(e),'/',int2str(MP.loops)]);
        combined{e} = zeros(tLength/binsize,MP.nStimuli);
        for s=1:nStim %MP.nStimuli % Bin spikes Alt: period = MP.transP_Test*MP.nTransPS; spikes(spikes >= (s-1)*period && spikes < s*period);
            spikes = trnESpikes{e,l}(sMembers{s},:).*MP.DT; % Input layer
            [stimSpikes{e,s},bins] = binSpikes(spikes,binsize,tLength);
            combined{e}(:,s) = stimSpikes{e,s}';
        end
    end
    
    if plotLevel > 0
    trainA = figure('Name', 'Training Input layer analysis');
    maximize(trainA);
    
    for e=1:lastEpoch
        disp(['[',simrun,'] ','Input layer training population analysis and autocorrelations, epoch ',int2str(e),'/',int2str(MP.loops)]);
        for s=1:nStim %MP.nStimuli
            subplot(nStim,3,((3*s)-2)); %MP.nStimuli
            plotFourierSpectrum(stimSpikes{e,s},binsize,trainA)
            title(['L',int2str(l-1),'S',int2str(s),' Fourier Spectrum']); % Power vs freq := periodogram
            
            subplot(nStim,3,(3*s)-1); %MP.nStimuli
            psth(stimSpikes{e,s},binsize,bins,tLength,MP.vExcit(l),trainA);
            title(['L',int2str(l-1),'S',int2str(s),' Excitatory Spike Histogram']);
            
            subplot(nStim,3,3*s) %MP.nStimuli
            nLags = MP.transP_Train / binsize; %MP.transP_Test / binsize; This was a bug?
            if useCorrFunc
                autocorr(stimSpikes{e,s},nLags); % Should this use binned spikes?
                set(gca, 'XTickMode', 'manual'); % Stops saveas recomputing XTicks
                set(gca,'XTickLabel',get(gca,'XTick')*binsize*1000); %Show ms not lags
            else
                [c_ww,lags]=xcorr(stimSpikes{e,s},nLags,'coeff');
                stem(lags(nLags+1:end)*binsize*1000,c_ww(nLags+1:end)); % +ve lags are enough
                ylabel('Correlation coefficient');
            end
            xlabel('Lag (ms)');
            title(['L',int2str(l-1),'S',int2str(s),' Autocorrelation']);
        end
        saveFig(trainA,['E',int2str(e-1),'L',int2str(l-1),'perStimHist',simrun],plotLevel);
        disp(['[',simrun,'] ','Done']);
        
        disp(['[',simrun,'] ','Plotting epoch ',int2str(m2c(e)),' Joint PSTH']);
        jPSTH=figure('Name', ['E',int2str(m2c(e)),'L',int2str(m2c(l)),' Joint PSTH']); %figure(jPSTH); 
        bar(bins*1000,combined{e}/(MP.vExcit(1)*binsize),1,'stacked'); % Modify psth?
        xlim([0,tLength*1000]);
        xlabel('Time (ms)');
        ylabel('Spikes/s');
        saveFig(jPSTH,['E',int2str(e-1),'L',int2str(l-1),'jointPSTH',simrun],plotLevel);
        close(jPSTH);
        disp(['[',simrun,'] ','Done']);
        
        %if e<MP.loops; close(trainA); end % Leave the last one
        
        if MP.nTestStimuli <= 4 % nchoosek is practical for n<~15
            disp(['[',simrun,'] ','Plotting epoch ',int2str(m2c(e)),' cross-correlations']);
            combs = nchoosek(1:MP.nTestStimuli,2);
            nCombs = size(combs,1);
            crossfh = figure('Name', ['E',int2str(m2c(e)),' cross-correlations']);
            maximize(crossfh);
            for c=1:nCombs
                subplot(nCombs,1,c);
                if useCorrFunc
                    crosscorr(stimSpikes{e,combs(c,1)},stimSpikes{e,combs(c,2)}); % Use last epoch: MP.loops not e
                    set(gca, 'XTickMode', 'manual'); % Stops saveas recomputing XTicks
                    if c==nCombs
                        xlabel('Lag (ms)');
                        set(gca,'XTickLabel',get(gca,'XTick')*binsize*1000); %Show ms not lags
                    else
                        xlabel('');
                        set(gca,'XTickLabel',[]);
                    end
                else
                    [c_ww,lags]=xcorr(stimSpikes{e,combs(c,1)},stimSpikes{e,combs(c,2)},nLags,'coeff');
                    stem(lags*binsize*1000,c_ww);
                end
                ylabel(['CCF: \{',int2str(combs(c,1)),',',int2str(combs(c,2)),'\}']);
                title('');
                %title(['E',int2str(e-1),'Inputs Cross-correlation: \{',int2str(combs(c,1)),',',int2str(combs(c,2)),'\}']);
            end
            saveFig(crossfh,['E',int2str(e-1),'inputsXcorr',simrun],plotLevel);
            close(crossfh);
            disp(['[',simrun,'] ','Done']);
        end
    end
    end
end

% Grun & Rotter 2010 Ch.2
% http://www.ton.scphys.kyoto-u.ac.jp/~shino/toolbox/english.htm

% optN = sshist(x); hist(x,optN);
% optW = sskernel(x); ksdensity(x,'width',optW);


%% Plot STDP curves

% This plots the shape of the STDP curves for Delta t ± 50ms according to 
% alphaC/D, tauC/D and learnR for a few initial weight values (Dg0). 
if plotLevel > 0
    disp(['[',simrun,'] ','Plotting STDP curves']);
    
    Dg0s = [0.25, 0.5, 0.75];
    t_post = 0;
    t_preN = 0.001 * (-50:-1);
    t_preP = 0.001 * (1:50);
    
    STDP = figure('Name', 'STDP curves');
    %maximize(STDP);
    labels = cell(length(Dg0s),1);
    hold all
    for w=1:length(Dg0s)
        Dg0 = Dg0s(w);
        labels{w} = ['\Deltag_0 = ',num2str(Dg0)];
        DgP = MP.learnR * (1-Dg0) * MP.alphaC*exp((t_preN-t_post)/MP.tauC);
        plot(1000*(t_post-t_preN),100*DgP/Dg0);%,'+-');
    end
    DgN = -MP.learnR * Dg0 * MP.alphaD*exp((t_post-t_preP)/MP.tauD);
    plot(1000*(t_post-t_preP),100*DgN/Dg0);%,'x-');
    % LTD is a constant '%' change (Dg0 cancels) - Perrinet et al., 2001 (Fig.1B)
    legend(labels); %,'Location','NorthWest')
    %plot(0,0,'+k','MarkerSize',10);
    plot(1000*[min(t_preN) max(t_preP)],[0,0],'k-');
    hold off
    xlabel('\Deltat = t_{post} - t_{pre} (msec)');
    ylabel('Percentage change in \Deltag');
    grid on; %minor;
    box on; 
    title(['Percentage weight change vs. time difference for ',int2str(length(Dg0s)),' initial weights']);
    % text('Position',[15,9],'String','Pre \rightarrow post: LTP', 'FontSize',14);
    % text('Position',[-45,-3],'String','Post \rightarrow pre: LTD', 'FontSize',14);
    saveFig(STDP,['STDPcurves',simrun],plotLevel);
    close(STDP);
    disp(['[',simrun,'] ','STDP curves plotted']);
end


%% Plot delays
if MP.axonDelay && MP.printConnections && plotLevel > 0
    
    ticks=4;
    barWidth=1; % Touching bars
    
    for l=1:MP.nLayers
        disp(['[',simrun,'] ','Layer ',int2str(l),' axonal delays']);
        
        if ~strcmp(MP.delayEfE,'MinD') && l>1
            dfh = figure('Name', ['L',int2str(m2c(l)),' EfE Delays']);
            
            % Plot Delay distribution
            subplot(1,2,1);
            delays = reshape(AffDelaysEfE{l-1}(:,2:end),1,[]);
            delays(delays==0) = []; % Remove padding zeros
            delays = delays*MP.DT*1000; % Convert to ms
            [freqs,dbins] = hist(delays,floor(sqrt(length(delays))));
            bar(dbins,freqs,'hist'); %barWidth
            %xlim([0 max(delays)]);
            xlabel('Axonal Delay bins (ms)');
            ylabel('Frequency');
            title({'Feed-Forward Axonal Delay Distribution'; ...
                ['Mean=',num2str(nanmean(delays),3),' SD=',num2str(nanstd(delays),3)]});
            
            % Plot Delay matrix
            subplot(1,2,2);
            delays = zeros(MP.vExcit(l),MP.vExcit(l-1));
            for n=1:MP.vExcit(l)
                nSyn = AffEfE{l-1}(n,1);
                delays(n,c2m(AffEfE{l-1}(n,1+(1:nSyn)))) = AffDelaysEfE{l-1}(n,1+(1:nSyn));
            end
            delays = delays*MP.DT*1000; % Convert to ms
            imagesc(delays);
            ch=colorbar;    ylabel(ch,'Delay (ms)');
            axis([1 MP.vExcit(l-1) 1 MP.vExcit(l)]);
            if (MP.vExcit(l-1) == MP.vExcit(l)); axis square; end
            set(gca,'YDir','normal'); % Normal y-axis
            set(gca,'YTick',MP.vExcit(l)/ticks:MP.vExcit(l)/ticks:MP.vExcit(l),...
                'XTick',MP.vExcit(l-1)/ticks:MP.vExcit(l-1)/ticks:MP.vExcit(l-1));
            %set(gca,'XTickMode','Manual'); %implicit
            xlabel('Input neuron');
            ylabel('Output neuron'); % Check these labels
            title('Feed-forward Axonal Delay Matrix');
            
            
            saveFig(dfh,['L',int2str(m2c(l)),'EfEDelays'],plotLevel);
            close(dfh);
        end
        
        if ~strcmp(MP.delayElE,'MinD') && (MP.pCnxElE(l) > eps)
            dfh = figure('Name', ['L',int2str(m2c(l)),' ElE Delays']);
            
            % Plot Delay distribution
            subplot(1,2,1);
            delays = reshape(AffDelaysElE{l}(:,2:end),1,[]);
            delays(delays==0) = []; % Remove padding zeros
            delays = delays*MP.DT*1000; % Convert to ms
            [freqs,dbins] = hist(delays,floor(sqrt(length(delays))));
            bar(dbins,freqs,'hist');
            %xlim([0 max(delays)]);
            xlabel('Axonal Delay bins');
            ylabel('Frequency');
            title({'Feed-Forward Axonal Delay Distribution'; ...
                ['Mean=',num2str(nanmean(delays),3),' SD=',num2str(nanstd(delays),3)]});
            
            % Plot Delay matrix
            subplot(1,2,2);
            delays = zeros(MP.vExcit(l),MP.vExcit(l));
            for n=1:MP.vExcit(l)
                nSyn = AffElE{l}(n,1);
                delays(n,c2m(AffElE{l}(n,1+(1:nSyn)))) = AffDelaysElE{l}(n,1+(1:nSyn));
            end
            delays = delays*MP.DT*1000; % Convert to ms
            imagesc(delays);
            ch=colorbar;    ylabel(ch,'Delay (ms)');
            axis([1 MP.vExcit(l) 1 MP.vExcit(l)]);
            axis square;
            set(gca,'YDir','normal'); % Normal y-axis
            set(gca,'YTick',MP.vExcit(l)/ticks:MP.vExcit(l)/ticks:MP.vExcit(l),...
                'XTick',MP.vExcit(l)/ticks:MP.vExcit(l)/ticks:MP.vExcit(l));
            %set(gca,'XTickMode','Manual'); %implicit
            xlabel('Input neuron');
            ylabel('Output neuron'); % Check these labels
            title('Lateral Axonal Delay Matrix');
            
            saveFig(dfh,['L',int2str(m2c(l)),'ElEDelays'],plotLevel);
            close(dfh);
        end
        
        if ~strcmp(MP.delayEI,'MinD') && (MP.pCnxEI(l) > eps)
            dfh = figure('Name', ['L',int2str(m2c(l)),' EI Delays']);
            
            % Plot Delay distribution
            subplot(1,2,1);
            delays = reshape(AffDelaysEI{l}(:,2:end),1,[]);
            delays(delays==0) = []; % Remove padding zeros
            delays = delays*MP.DT*1000; % Convert to ms
            [freqs,dbins] = hist(delays,floor(sqrt(length(delays))));
            bar(dbins,freqs,'hist');
            %xlim([0 max(delays)]);
            xlabel('Axonal Delay bins');
            ylabel('Frequency');
            title({'EI Axonal Delay Distribution'; ...
                ['Mean=',num2str(nanmean(delays),3),' SD=',num2str(nanstd(delays),3)]});
            
            % Plot Delay matrix
            subplot(1,2,2);
            delays = zeros(MP.vInhib(l),MP.vExcit(l));
            for n=1:MP.vInhib(l)
                nSyn = AffEI{l}(n,1);
                delays(n,c2m(AffEI{l}(n,1+(1:nSyn)))) = AffDelaysEI{l}(n,1+(1:nSyn));
            end
            delays = delays*MP.DT*1000; % Convert to ms
            imagesc(delays);
            ch=colorbar;    ylabel(ch,'Delay (ms)');
            axis([1 MP.vExcit(l) 1 MP.vInhib(l)]);
            if (MP.vExcit(l) == MP.vInhib(l)); axis square; end
            set(gca,'YDir','normal'); % Normal y-axis
            set(gca,'YTick',MP.vInhib(l)/ticks:MP.vInhib(l)/ticks:MP.vInhib(l),...
                'XTick',MP.vExcit(l)/ticks:MP.vExcit(l)/ticks:MP.vExcit(l));
            %set(gca,'XTickMode','Manual'); %implicit
            xlabel('Input neuron');
            ylabel('Output neuron'); % Check these labels
            title('EI Axonal Delay Matrix');
            
            saveFig(dfh,['L',int2str(m2c(l)),'EIDelays'],plotLevel);
            close(dfh);
        end
        disp(['[',simrun,'] ','Layer ',int2str(l),' Done!']);
    end
end


%% Scatter plot and distribution of lateral weights

ticks=4;
% if plotLevel && (MP.trainElE || MP.SOM) % Adapt for plotting Gaussian
% lateral weight profile
% Also adapt to plot lateral delays (and FF delays)
if MP.trainElE && plotLevel > 0
    
    if MP.priorPhases % Reconsider.........................................
        LatWeights = PP_WeightsElE;
        if MP.pretrain
            ptLatWeights = ptPP_WeightsElE;
        end
    else
        LatWeights = WeightsElE;
        if MP.pretrain
            ptLatWeights = ptWeightsElE;
        end
    end
    
    for l=1:MP.nLayers
        if (MP.pCnxElE(l) > eps)
            if MP.pretrain
                ptlatWeights = figure('Name', ['Pre Training L',int2str(m2c(l)),' Lateral Weights']);
                subplot(2,2,1);
                c=0;
                x=zeros(1,MP.vExcit(l)^2);
                y=zeros(1,MP.vExcit(l)^2);
                %cx=1; cy=1;
                for n=1:MP.vExcit(l)
                    nSyn = AffElE{l}(n,1);
                    vSyn = c2m(AffElE{l}(n,1+(1:nSyn)));
                    for s=1:nSyn
                        preN = vSyn(s);
                        preS = find(AffElE{l}(preN,2:end)==m2c(n));
                        if ~isempty(preS) % Extract Bi-directional weights
                            c = c + 1;
                            x(c) = ptLatWeights{l}(n,s);
                            y(c) = ptLatWeights{l}(preN,preS);
                        else % Extract uni-directional weight
                            %WeightsElE{l}(n,s);
                        end
                    end
                end
                scatter(x(1:c),y(1:c));
                hold on;
                plot([0,1],[0,1],'b:');
                hold off;
                axis([0 1 0 1]); axis square; box on;
                xlabel('$\Delta g_{ij}$','interpreter','latex');
                ylabel('$\Delta g_{ji}$','interpreter','latex');
                title(['L',int2str(m2c(l)),' bi-directional lateral weights']);
                
                subplot(2,2,3)
                weights = zeros(MP.vExcit(l),MP.vExcit(l));
                for n=1:MP.vExcit(l)
                    nSyn = AffElE{l}(n,1);
                    weights(n,c2m(AffElE{l}(n,1+(1:nSyn)))) = ptLatWeights{l}(n,1:nSyn);
                end
                if MP.stimGroups && l==1
                    sVec = zeros(1,MP.vExcit(l));
                    offset = 0;
                    for g=1:MP.nGroups
                        nonshared = gMembers{g}(~ismember(gMembers{g},gShared));
                        sVec((1:length(nonshared))+offset) = nonshared;
                        offset = offset + length(nonshared);
                    end
                    sVec((1:length(gShared))+offset) = gShared;
                    weights = weights(sVec,:);
                    weights = weights(:,sVec);
                end
                imagesc(weights);%WeightsElE{l});
                set(gca,'YDir','normal'); % Normal y-axis
                set(gca,'YTick',MP.vExcit(l)/ticks:MP.vExcit(l)/ticks:MP.vExcit(l),...
                    'XTick',MP.vExcit(l)/ticks:MP.vExcit(l)/ticks:MP.vExcit(l));
                axis([1 MP.vExcit(l) 1 MP.vExcit(l)]);
                axis square;
                xlabel('Input neuron');
                ylabel('Output neuron'); % Check these labels
                title('Synaptic Weight Matrix');
                
                subplot(2,2,[2 4])
                weights = reshape(ptLatWeights{l},1,[]);
                %weights(isnan(weights))=[];
                hist(weights,floor(sqrt(length(weights))));
                xlim([0.0 1.0]);
                xlabel('Synaptic weight bins');
                ylabel('Frequency');
                title({'Synaptic Weight Distribution'; ...
                    ['Mean=',num2str(nanmean(weights),3),' SD=',num2str(nanstd(weights),3)]});
                
                saveFig(ptlatWeights,['ptL',int2str(m2c(l)),'LateralWeights'],plotLevel);
            end
            
            
            latWeights = figure('Name', ['Post Training L',int2str(m2c(l)),' Lateral Weights']);
            subplot(2,2,1);
            c=0;
            x=zeros(1,MP.vExcit(l)^2);
            y=zeros(1,MP.vExcit(l)^2);
            for n=1:MP.vExcit(l)
                nSyn = AffElE{l}(n,1);
                vSyn = c2m(AffElE{l}(n,1+(1:nSyn)));
                for s=1:nSyn
                    preN = vSyn(s);
                    preS = find(AffElE{l}(preN,2:end)==m2c(n));
                    if ~isempty(preS) % Extract Bi-directional weights
                        c = c + 1;
                        x(c) = LatWeights{l}(n,s);
                        y(c) = LatWeights{l}(preN,preS);
                    else % Extract uni-directional weight
                        %WeightsElE{l}(n,s);
                    end
                end
            end
            scatter(x(1:c),y(1:c));
            hold on;
            plot([0,1],[0,1],'b:');
            hold off;
            axis([0 1 0 1]); axis square; box on;
            xlabel('$\Delta g_{ij}$','interpreter','latex');
            ylabel('$\Delta g_{ji}$','interpreter','latex');
            title(['L',int2str(m2c(l)),' bi-directional lateral weights']);
            
            subplot(2,2,3)
            weights = zeros(MP.vExcit(l),MP.vExcit(l));
            for n=1:MP.vExcit(l)
                nSyn = AffElE{l}(n,1);
                weights(n,c2m(AffElE{l}(n,1+(1:nSyn)))) = LatWeights{l}(n,1:nSyn);
            end
            if MP.stimGroups && l==1
                sVec = zeros(1,MP.vExcit(l));
                offset = 0;
                for g=1:MP.nGroups
                    nonshared = gMembers{g}(~ismember(gMembers{g},gShared));
                    sVec((1:length(nonshared))+offset) = nonshared;
                    offset = offset + length(nonshared);
                end
                sVec((1:length(gShared))+offset) = gShared;
                weights = weights(sVec,:);
                weights = weights(:,sVec);
            end
            imagesc(weights);%WeightsElE{l});
            set(gca,'YDir','normal'); % Normal y-axis
            set(gca,'YTick',MP.vExcit(l)/ticks:MP.vExcit(l)/ticks:MP.vExcit(l),...
                'XTick',MP.vExcit(l)/ticks:MP.vExcit(l)/ticks:MP.vExcit(l));
            axis([1 MP.vExcit(l) 1 MP.vExcit(l)]);
            axis square;
            xlabel('Input neuron');
            ylabel('Output neuron'); % Check these labels
            title('Synaptic Weight Matrix');
            
            subplot(2,2,[2 4])
            weights = reshape(LatWeights{l},1,[]);
            %weights(isnan(weights))=[]; %Unnecessary
            hist(weights,floor(sqrt(length(weights))));
            xlim([0.0 1.0]);
            xlabel('Synaptic weight bins');
            ylabel('Frequency');
            title({'Synaptic Weight Distribution'; ...
                ['Mean=',num2str(nanmean(weights),3),' SD=',num2str(nanstd(weights),3)]});
            
            saveFig(latWeights,['L',int2str(m2c(l)),'LateralWeights'],plotLevel);
        end
    end
end


%% Plot Rasters and weights

if plotLevel > 0
% Calculate weight bins and centres
%wbins = 0.025:0.05:0.975;
%gCols = colormap(hsv(MP.nGroups));
linecols=colormap(lines);
colormap('default'); % Default is jet(64)
%tmp = colormap;
%gCols = tmp(linspace(1,size(tmp,1),MP.nGroups),:);
if MP.stimGroups
    nCols = MP.nGroups;
    if ~isempty(gShared) %gShared = find(sum(prototypes)>1);
        nCols = nCols + 1;
    end
    gCols = linecols(1:nCols,:);
end


if MP.pretrain==1
    
    % Plot pretraining inputs 
    if (plotInput && MP.trainElE)
        
        if MP.priorPhases
            disp(['[',simrun,'] ','Plotting post prior training inputs']);
            postPriorInputs = figure('Name', 'Plotting post prior training inputs');
            PPtestMS = MP.nStimuli * MP.nTransPS * MP.transP_Test; %EpochMS; % Caution!
            spikes = PP_ESpikes{1}*MP.DT*1000;
            cSpec = zeros(size(spikes,1),3); % Colour specified as [r g b] triplets
            for g=1:MP.nGroups
                members = find(prototypes(g,:));
                cSpec(members,:) = repmat(gCols(g,:),[length(members) 1]);
            end
            cSpec(gShared,:) = repmat(gCols(nCols,:),[length(gShared) 1]);
            raster(spikes,PPtestMS,true,cSpec);
            hold on; plot(zeros(size(gShared)),gShared,'*','Color',gCols(nCols,:)); hold off;
            if ~MP.useFilteredImages %isfield(MP,'nFiringNeurons')
                set(gca,'YTick',0:MP.nFiringNeurons:MP.vExcit(1));
            end
            if MP.nTransPS>1; set(gca,'XTick',0:(1000*MP.transP_Test):PPtestMS); end
            saveFig(postPriorInputs,['postPriorInputs',simrun],plotLevel);
            close(postPriorInputs);
            disp(['[',simrun,'] ','Done']);
        end
        
        disp(['[',simrun,'] ','Plotting pre-training inputs']);
        preInputs = figure('Name', 'Plotting pre-training inputs');
        spikes = ptESpikes{1}*MP.DT*1000;
        if (MP.stimGroups)
            cSpec = zeros(size(spikes,1),3); % Colour specified as [r g b] triplets
            for g=1:MP.nGroups
                members = find(prototypes(g,:));
                cSpec(members,:) = repmat(gCols(g,:),[length(members) 1]);
            end
            cSpec(gShared,:) = repmat(gCols(nCols,:),[length(gShared) 1]);
            raster(spikes,TestMS,true,cSpec);
            hold on; plot(zeros(size(gShared)),gShared,'*','Color',gCols(nCols,:)); hold off;
        else
            raster(spikes,TestMS);%,MP.DT); %L0ExcitSpikes %ptESpikes{1}*MP.DT*1000
        end
        %set(gca,'YTick',0:MP.vExcit(MP.nLayers)/8:MP.vExcit(MP.nLayers));
        if ~MP.useFilteredImages %isfield(MP,'nFiringNeurons')
            set(gca,'YTick',0:MP.nFiringNeurons:MP.vExcit(1));
        end
        if MP.nTestTransPS>1; set(gca,'XTick',0:(1000*MP.transP_Test):TestMS); end
        saveFig(preInputs,['ptInputs',simrun],plotLevel);
        close(preInputs);
        disp(['[',simrun,'] ','Done']);
    end
    
    % Plot pre-training results
    disp(['[',simrun,'] ','Plotting pretraining output raster and weight distribution']);
    preAnalysis = figure('Name', 'Pretraining output raster and weight distribution');
    if MP.nLayers>1
        subplot(2,3,1);
        if MP.printConnections
            weights = zeros(MP.vExcit(MP.nLayers),MP.vExcit(MP.nLayers-1));
            for n=1:MP.vExcit(MP.nLayers)
                nSyn = AffEfE{MP.nLayers-1}(n,1);
                weights(n,c2m(AffEfE{MP.nLayers-1}(n,1+(1:nSyn)))) = ptWeightsEfE{MP.nLayers-1}(n,1:nSyn);
            end
        else % Legacy
            weights = ptWeightsEfE{MP.nLayers-1}; 
        end
        imagesc(weights);
        axis([1 MP.vExcit(MP.nLayers-1) 1 MP.vExcit(MP.nLayers)]);
        if (MP.vExcit(MP.nLayers-1) == MP.vExcit(MP.nLayers)); axis square; end
        set(gca,'YDir','normal'); % Normal y-axis
        set(gca,'YTick',MP.vExcit(MP.nLayers)/ticks:MP.vExcit(MP.nLayers)/ticks:MP.vExcit(MP.nLayers),...
            'XTick',MP.vExcit(MP.nLayers-1)/ticks:MP.vExcit(MP.nLayers-1)/ticks:MP.vExcit(MP.nLayers-1));
        %set(gca,'XTickMode','Manual'); %Implicit
        xlabel('Input neuron');
        ylabel('Output neuron'); % Check these labels
        title('Synaptic Weight Matrix');
        
        subplot(2,3,4);
        weights = reshape(ptWeightsEfE{MP.nLayers-1},1,[]);
        hist(weights,floor(sqrt(length(weights)))); % hist(reshape(ptWeights{MP.nLayers-1},1,[]),wbins);
        xlim([0 1]);
        xlabel('Synaptic weight bins');
        ylabel('Frequency');
        title({'Synaptic Weight Distribution'; ...
            ['Mean=',num2str(nanmean(weights),3),' SD=',num2str(nanstd(weights),3)]});
        
        subplot(2,3,[2 6]);
    end
    
    spikes = ptESpikes{MP.nLayers}*MP.DT*1000;
    if (MP.stimGroups && MP.nLayers == 1)
        cSpec = zeros(size(spikes,1),3); % Colour specified as [r g b] triplets
        for g=1:MP.nGroups
            members = find(prototypes(g,:));
            cSpec(members,:) = repmat(gCols(g,:),[length(members) 1]);
        end
        cSpec(gShared,:) = repmat(gCols(nCols,:),[length(gShared) 1]);
        raster(spikes,TestMS,true,cSpec);
        hold on; plot(zeros(size(gShared)),gShared,'*','Color',gCols(nCols,:)); hold off;
    else
        raster(spikes,TestMS);%,MP.DT); %pt_output_spikes %ptESpikes{MP.nLayers}*MP.DT*1000
    end
    set(gca,'YTick',0:MP.vExcit(MP.nLayers)/8:MP.vExcit(MP.nLayers));
    if MP.nTestTransPS>1; set(gca,'XTick',0:(1000*MP.transP_Test):TestMS); end
    saveFig(preAnalysis,['ptOutputAnalysis',simrun],plotLevel);
    disp(['[',simrun,'] ','Done']);
end


% Plot training stimuli (first layer only)
if plotTraining==true && MP.train==true
    if (plotFullTraining)
        lastEpoch=MP.loops-1;
    else
        lastEpoch=0;
    end
    for e=0:lastEpoch
        disp(['[',simrun,'] ','Plotting epoch ',int2str(e),' training stimuli']);
        trainInputs = figure('Name', ['E',int2str(e),' training stimuli']);
        maximize(trainInputs);
        spikes = trnESpikes{e+1,1}*MP.DT*1000;
        if (MP.stimGroups)
            cSpec = zeros(size(spikes,1),3); % Colour specified as [r g b] triplets
            for g=1:MP.nGroups
                members = find(prototypes(g,:));
                cSpec(members,:) = repmat(gCols(g,:),[length(members) 1]);
            end
            cSpec(gShared,:) = repmat(gCols(nCols,:),[length(gShared) 1]);
            raster(spikes,EpochMS,true,cSpec);
            hold on; plot(zeros(size(gShared)),gShared,'*','Color',gCols(nCols,:)); hold off;
        else
            raster(spikes,EpochMS);%,MP.DT); % Only plotting inputs %trnESpikes{e+1,1}*MP.DT*1000
        end
        if ~MP.useFilteredImages %isfield(MP,'nFiringNeurons')
            set(gca,'YTick',0:MP.nFiringNeurons:MP.vExcit(1));
        end
        if MP.nTransPS>1; set(gca,'XTick',0:(1000*MP.transP_Train):EpochMS); end
        grid('minor');
        saveFig(trainInputs,['E',int2str(e),'inputs',simrun],plotLevel);
        if e<MP.loops-1; close(trainInputs); end % Leave the last one
        disp(['[',simrun,'] ','Done']);
    end
end

% Plot post-training results
if plotInput==true
    disp(['[',simrun,'] ','Plotting post-training inputs']);
    postInputs = figure('Name', 'Post-training inputs');
    spikes = ESpikes{1}*MP.DT*1000;
    if (MP.stimGroups)
        cSpec = zeros(size(spikes,1),3); % Colour specified as [r g b] triplets
        for g=1:MP.nGroups
            %set(gcf,'Color',idmap(g,:));
            %set(gcf,'DefaultAxesColorOrder',[1 0 0;0 0 1;0 1 0])
            members = find(prototypes(g,:));
            cSpec(members,:) = repmat(gCols(g,:),[length(members) 1]);
            %spec='-b';
        end
        cSpec(gShared,:) = repmat(gCols(nCols,:),[length(gShared) 1]);
        raster(spikes,TestMS,true,cSpec);
        hold on; plot(zeros(size(gShared)),gShared,'*','Color',gCols(nCols,:)); hold off;
    else
        raster(spikes,TestMS);%,MP.DT); %L0ExcitSpikes %ESpikes{1}*MP.DT*1000
    end
    %set(gca,'YTick',0:MP.vExcit(MP.nLayers)/8:MP.vExcit(MP.nLayers));
    if ~MP.useFilteredImages %isfield(MP,'nFiringNeurons')
        set(gca,'YTick',0:MP.nFiringNeurons:MP.vExcit(1));
    end
    if MP.nTestTransPS>1; set(gca,'XTick',0:(1000*MP.transP_Test):TestMS); end
    saveFig(postInputs,['inputs',simrun],plotLevel);
    close(postInputs);
    disp(['[',simrun,'] ','Done']);
end

if plotInhib == 1
    disp(['[',simrun,'] ','Plotting Inhibitory spike rasters']);
    if MP.inputInhib == 1
        inpInhib = figure('Name', 'Inhibitory Input layer spike rasters');
        raster(ISpikes{1}*MP.DT*1000,TestMS);%,MP.DT); %L0InhibSpikes
        saveFig(inpInhib,['InhibInputs',simrun],plotLevel);
        close(inpInhib);
    end

    % Plot output layer inhibitory spikes
    outInhib = figure('Name', 'Inhibitory Output layer spike rasters');
    raster(ISpikes{MP.nLayers}*MP.DT*1000,TestMS);%,MP.DT); %eval(['L',int2str(nLayers-1),'InhibSpikes'])
    saveFig(outInhib,['InhibOutputs',simrun],plotLevel);
    disp(['[',simrun,'] ','Done']);
end

disp(['[',simrun,'] ','Plotting post-training output raster and weight distribution']);
postAnalysis = figure('Name', 'Post-training output raster and weight distribution');
if MP.nLayers>1
subplot(2,3,1);
if MP.printConnections
    weights = zeros(MP.vExcit(MP.nLayers),MP.vExcit(MP.nLayers-1));
    for n=1:MP.vExcit(MP.nLayers)
        nSyn = AffEfE{MP.nLayers-1}(n,1);
        weights(n,c2m(AffEfE{MP.nLayers-1}(n,1+(1:nSyn)))) = WeightsEfE{MP.nLayers-1}(n,1:nSyn);
    end
else % Legacy
    weights = WeightsEfE{MP.nLayers-1};
end
imagesc(weights);
axis([1 MP.vExcit(MP.nLayers-1) 1 MP.vExcit(MP.nLayers)]);
if (MP.vExcit(MP.nLayers-1) == MP.vExcit(MP.nLayers)); axis square; end
set(gca,'YDir','normal'); % Normal y-axis
set(gca,'YTick',MP.vExcit(MP.nLayers)/ticks:MP.vExcit(MP.nLayers)/ticks:MP.vExcit(MP.nLayers),...
        'XTick',MP.vExcit(MP.nLayers-1)/ticks:MP.vExcit(MP.nLayers-1)/ticks:MP.vExcit(MP.nLayers-1));
%set(gca,'XTickMode','Manual'); %implicit
xlabel('Input neuron');
ylabel('Output neuron'); % Check these labels
title('Synaptic Weight Matrix');

subplot(2,3,4);
weights = reshape(WeightsEfE{MP.nLayers-1},1,[]);
hist(weights,floor(sqrt(length(weights)))); %hist(reshape(Weights{MP.nLayers-1},1,[]),wbins);
xlim([0 1]);
xlabel('Synaptic weight bins');
ylabel('Frequency');
title({'Synaptic Weight Distribution'; ...
    ['Mean=',num2str(nanmean(weights),3),' SD=',num2str(nanstd(weights),3)]});

subplot(2,3,[2 6]);
end
raster(ESpikes{MP.nLayers}*MP.DT*1000,TestMS);%,MP.DT); %eval(['L',int2str(nLayers-1),'ExcitSpikes'])
set(gca,'YTick',0:MP.vExcit(MP.nLayers)/8:MP.vExcit(MP.nLayers));
if MP.nTransPS>1; set(gca,'XTick',0:(1000*MP.transP_Train):TestMS); end
saveFig(postAnalysis,['OutputAnalysis',simrun],plotLevel);
disp(['[',simrun,'] ','Done']);

end


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
%clear remain token

if (savePDF)
    system('mkdir eps');
    system('mv *.eps eps');
    system('mkdir pdf');
    system('mv *.pdf pdf');
end

% if plotLevel && ~(MP.nRecords && plotRecords)
%     %close all;
%     system('open *.png'); % Opens all pngs in the same window in 10.7
% end

%% Cell Analysis

% regime = 'Testing';%'Training';
% if (strcmpi(regime, 'Training')) && MP.loops > 1
%     loops=[0,MP.loops-1]; % Plot first and last
% else
%     loops = 0;
% end

pCount = 1; % Always a testing phase
% if MP.priorPhases
%     if MP.pretrain
%         pCount = pCount + 2;
%     end
%     if MP.train
%         pCount = pCount + 1;
%     end
% end
if MP.pretrain
    pCount = pCount + 1;
end
if MP.train
    pCount = pCount + 1;
end

phases = cell(pCount,1);

%if MP.priorPhases... % PP data will not work without modifying code below

pCount = 0;
if MP.pretrain
    pCount = pCount + 1;
    phases{pCount} = 'preTraining';
end
if MP.train
    pCount = pCount + 1;
    phases{pCount} = 'Training';
end
pCount = pCount + 1;
phases{pCount} = 'Testing';

%ticks = 0:0.25:1;
    
if (MP.nRecords && plotRecords) %MP.nRecordsPL
    disp(['[',simrun,'] ','Plotting individual cell records']);
    records;
    
    %if exist('MP.Records', 'var')
    %    d = '';
    if exist('records.tbz','file')
        system('tar -xvf records.tbz'); 
        d = '';
    elseif exist('Records','dir') % Legacy
        d = 'Records/';
    end
%    system('mv R*.dat Records');

if incInj
injTest = zeros(MP.nTestStimuli, MP.nTestStimuli * MP.nTestTransPS);
for s=1:MP.nTestStimuli
    injTest(s,(1:MP.nTestTransPS)+(s-1)*MP.nTestTransPS) = MP.current;
end

if MP.train
    injTrain = cell(MP.loops,1);
    for e=1:MP.loops
        injTrain{e} = zeros(MP.nStimuli, MP.nStimuli*MP.nTransPS);
        for s=1:MP.nStimuli % Use schedule to calc when each transform was shown
            injTrain{e}(s,(STIM.schedule{e}(1,:)==m2c(s))) = MP.current;
        end
    end
end

end

for p = 1:pCount
    if (strcmpi(phases{p}, 'preTraining'))
        loops = 0;
        timeMS = TestMS;
        weights = ptWeightsEfE; % Pre training
        nStim = MP.nTestStimuli;
        nTrans = MP.nTestTransPS;
        %if MP.newTestSet; pats = STIM.test; else pats = STIM.train; end
        if incInj %exist('STIM','var')
            pats = STIM.test;
            transMS = 1000*MP.transP_Test;
            sInj = zeros(nStim, 1+timeMS);
            for s=1:nStim
                sInj(s,2:end) = reshape(repmat(injTest(s,:),[transMS,1]),1,[]);
            end
        end
        
    elseif (strcmpi(phases{p}, 'Training'))
        if MP.loops > 1
            loops=[0,MP.loops-1]; % Plot first and last
        else
            loops = 0;
        end
        timeMS = EpochMS; %RecordsMS;
        % Plot evolution of weights video.
        weights = WeightsEfE; % Tempory placeholder...
        nStim = MP.nStimuli;
        nTrans = MP.nTransPS;
        if incInj
            pats = STIM.train;
            transMS = 1000*MP.transP_Train;
        end

    elseif (strcmpi(phases{p}, 'Testing'))
        loops = 0;
        timeMS = TestMS;
        weights = WeightsEfE; % Post training
        nStim = MP.nTestStimuli;
        nTrans = MP.nTestTransPS;
        %if MP.newTestSet; pats = STIM.test; else pats = STIM.train; end
        if incInj
            pats = STIM.test;
            transMS = 1000*MP.transP_Test;
            sInj = zeros(nStim, 1+timeMS);
            for s=1:nStim
                sInj(s,2:end) = reshape(repmat(injTest(s,:),[transMS,1]),1,[]);
            end
        end
        
    else
        disp('Error: Unknown regime!');
    end
    for l=0:MP.nLayers-1
        for r=1:MP.vRecords(c2m(l))%MP.nRecordsPL
            nInj = [];
            n = MP.Records{c2m(l)}(r);
            for e=1:length(loops)
                if (strcmpi(phases{p}, 'preTraining'))
                    spikes = ptESpikes{c2m(l)}(c2m(n),:)*MP.DT*1000;
                    prefix = [d,'RptL',int2str(l),'N',int2str(n)]; %Records/
                    if incInj && l==0
                        stMem = find(STIM.protoTest(:,c2m(n)));
                        nInj = zeros(1, nStim * nTrans);
                        %for s=1:length(stMem)
                        for t=1:nTrans
                            nInj(t+(stMem-1)*nTrans) = pats{stMem,t}(c2m(n)); %stMem(s)
%                             if STIM.test{stMem,t}(c2m(n)) > 0
%                                 nInj(t+(stMem-1)*MP.nTestTransPS) = MP.current;
%                             end
                        end
                        %end
                        nInj = reshape(repmat(nInj,[transMS,1]),1,[]);
                        nInj = [0,nInj];
                    end
                    %else
%                         injTest = zeros(MP.nTestStimuli, MP.nTestStimuli * MP.nTestTransPS);
%                         for s=1:MP.nTestStimuli
%                             injTest(s,(1:MP.nTestTransPS)+(s-1)*MP.nTestTransPS) = MP.current;
%                         end
                    %end
%                     sInj = zeros(nStim, 1+timeMS);
%                     for s=1:MP.nStim
%                         sInj(s,2:end) = reshape(repmat(injTest(s,:),[transMS,1]),1,[]);
%                     end
                        
                elseif (strcmpi(phases{p}, 'Training'))
                    spikes = trnESpikes{c2m(loops(e)),c2m(l)}(c2m(n),:)*MP.DT*1000;
                    prefix = [d,'RE',int2str(loops(e)),'L',int2str(l),'N',int2str(n)];
                    if incInj % Pass injected current
                        if l==0 % Schedule for specific neuron
%                             stMem = find(STIM.protoTrain(:,c2m(n)));
                            nInj = zeros(1, nStim * nTrans);
                            pat=0;
                            for s=1:nStim
                                for t=1:nTrans
                                    pat=pat+1;
                                    st = c2m(STIM.schedule{c2m(loops(e))}(1,pat));
                                    tr = c2m(STIM.schedule{c2m(loops(e))}(2,pat));
                                    nInj(pat) = pats{st,tr}(c2m(n));
                                end
                            end
                            nInj = reshape(repmat(nInj,[transMS,1]),1,[]);
                            nInj = [0,nInj];
                        end
%                         else % plot each stimulus's current
%                             injTrain = zeros(MP.nStimuli, MP.nStimuli*MP.nTransPS);
%                             for s=1:MP.nStimuli % Use schedule to calc when each transform was shown
%                                 injTrain(s,(STIM.schedule{c2m(l)}(1,:)==m2c(s))) = MP.current;
%                             end
%                         end
                        
                        sInj = zeros(nStim, 1+timeMS);
                        for s=1:nStim % Reconsider for K>1
                            sInj(s,2:end) = reshape(repmat(injTrain{e}(s,:),[transMS,1]),1,[]);
                        end
                    end
                    
                elseif (strcmpi(phases{p}, 'Testing'))
                    spikes = ESpikes{c2m(l)}(c2m(n),:)*MP.DT*1000;
                    prefix = [d,'RL',int2str(l),'N',int2str(n)];
                    if incInj && l==0
                        stMem = find(STIM.protoTest(:,c2m(n)));
                        nInj = zeros(1, nStim * nTrans);
                        %for s=1:length(stMem)
                        for t=1:nTrans
                            nInj(t+(stMem-1)*nTrans) = pats{stMem,t}(c2m(n));
                        end
                        %end
                        nInj = reshape(repmat(nInj,[transMS,1]),1,[]);
                        nInj = [0,nInj];
                    end
                else
                    disp('Error: Unknown regime!');
                end
                
                if incInj
                    ca = cellAnalysis(MP,l,n,spikes,prefix,loops(e),timeMS,sInj,nInj);
                else
                    ca = cellAnalysis(MP,l,n,spikes,prefix,loops(e),timeMS);
                end
                maximize(ca);
                set(gcf,'PaperPositionMode', 'manual', ...
                    'PaperType', 'A4', 'PaperUnits','centimeters', ...
                    'Paperposition',[1 1 28.7 20]);
                orient(ca,'landscape'); % 'PaperOrientation', 'landscape',
                saveFig(ca,[prefix,'cellAnalysis',simrun],plotLevel); % *** Save as A4 landscape
                %saveas(ca, [prefix,'cellAnalysis',simrun],'png');
                close(ca);
            end
            if plotArchitecture && l>0 && MP.printConnections % Efferent connections not yet implemented
            % 3D plot - Pass figure handle and plot all connections to/from a particular neuron
            cnx3D = figure('Name', 'Network structure');
            plotNet3D(MP, weights, AffEfE, cnx3D, l, n); % Check this is neuron ID not record ID
            title('Excitatory feed-forward post-training weights'); %---<<<
            %maximize(cnx3D);
            % Change camera angle
            saveFig(cnx3D, [prefix,'connectivity',simrun],plotLevel);
            %saveas(cnx3D, [prefix,'connectivity',simrun],'png'); % Post training
            close(cnx3D);
            end
        end
    end
end
    disp(['[',simrun,'] ','Cell records plotted']);
%     if plotLevel
%         %close all;
%         system(['open ',d,'R*.png']); % Opens all pngs in the same window in 10.7
%     end
    if (plotArchitecture && MP.printConnections)
        disp(['[',simrun,'] ','Plotting whole network structure']);
        fullNet = plotNet3D(MP, weights, AffEfE); % Full network plot
        %maximize(fullNet);
        title('Network Connectivity');
        %saveas(fullNet, ['NetConnectivity',simrun],'png'); % Post training
        disp(['[',simrun,'] ','Done']);
    end
    %if exist('records.tbz','file'); system('rm *.dat'); end
    if exist('records.tbz','file'); delete([d,'R*.dat']); end
end

if plotLevel
    system('open *.png'); % Opens all pngs in the same window in 10.7
end

%% Clean up
clear remain token Dg0 SOMfRates STDP ans e fRates fRdist fname freqA l loop
clear multInfo singInfo n postAnalysis postInputs prefix st st2 tr
clear t_post t_preN t_preP DgN DgP trainInputs

save matlab_workspace
disp(['[',simrun,'] ','Results saved in ', '<a href = "file://',pwd,filesep,'matlab_workspace.mat">',strcat(rdir,simrun),filesep,'matlab_workspace.mat</a>']);
toc;


% function index = c2m(index) % C --> Matlab indexing
% index=index+1;


% function index = m2c(index) % Matlab --> C indexing
% index=index-1;

% Extend these to deal with row major (C) <--> column major (Matlab)?


% function [freqs,bins] = binSpikes(spikes, binsize, tLength)
% %spikes(spikes(:)==0)=[];
% spikes=spikes(:)'; % Converts matrix to vector % reshape(A, 1, prod(size(A)))
% %spikes(isnan(spikes(:)))=[];
% spikes(spikes(:)<0)=[]; % Not strictly necessary since bins start at 0
% bins = 0:binsize:tLength;
% freqs = histc(spikes,bins);
% freqs(end) = [];
% bins(end) = [];
% %freqs = freqs/(nTrials*binsize); % Normalise ???


% function plotFourierSpectrum(freqs,binsize,fh)%tLength,fh)
% % spikes: vector of spike times in seconds
% % binsize: Bin size in seconds
% % tLength: Total recorded time
% % fh: Figure handle
% 
% % http://www.mathworks.com/products/matlab/demos.html?file=/products/demos/shipping/matlab/sunspots.html
% 
% % bins = binsize/2 : binsize : (tLength - binsize/2); % specify centres for hist(), or edges for histc()
% % freqs = hist(spikes,bins);
% 
% %freqs = freqs - mean(freqs);
% 
% % bins = 0:binsize:tLength;
% % freqs = histc(spikes,bins);
% % freqs(end) = [];
% 
% Fs = 1/binsize; % Sampling frequency (Hz)
% NyquistLim = Fs/2; % (1/sampInt)/2
% disp(['Nyquist Limit = ',num2str(NyquistLim,3),' Hz']);
% 
% NFFT = 2^nextpow2(length(freqs)); % Next power of 2 from length of freqs %http://blinkdagger.com/matlab/matlab-fft-and-zero-padding/
% FFT = fft(freqs,NFFT);%/length(freqs);%Fs;
% FFT(1) = []; % This is simply the sum of the data and can be removed
% % Plot first half of FFT (second half is a reflection of the first)
% f = NyquistLim * linspace(0,1,NFFT/2 + 1); % Gives f up to Nyquist limit NyquistLim * (0:NFFT-1); %
% %figure(fh); plot(f,2*abs(FFT(1:NFFT/2 +1))/NFFT);    ylabel('|Amplitute|');
% figure(fh); plot(f,(abs(FFT(1:NFFT/2 + 1)).^2)/NFFT);   ylabel('Power'); % Periodogram
% xlabel('Frequency (Hz)');
% xlim([0,NyquistLim]); %xlim([-1,NyquistLim+1]);


% function psth(freqs,binsize,bins,tLength,nTrials,fh)
% % spikes: vector of spike times in seconds
% % binsize: Bin size in seconds
% % tLength: Total recorded time
% % nTrials: The number of trials (or cells) to average over
% % fh: Figure handle
% 
% % bins = binsize/2 : binsize : (tLength - binsize/2); % specify centres for hist(), or edges for histc()
% % freqs = hist(spikes,bins)/(nTrials*binsize); % Normalise
% %freqs = freqs - mean(freqs);
% freqs = freqs/(nTrials*binsize); % Normalise
% figure(fh); bar(bins*1000,freqs,1); % Plot as msec
% xlim([0,tLength*1000]); %+(step/2)
% xlabel('Time (ms)'); %Seconds (Binned)
% %ylabel('Average Spikes/s'); % Without normalising would just be frequency
% ylabel('Spikes/s');
% 
% % Grun & Rotter 2010 Ch.2
% % http://www.ton.scphys.kyoto-u.ac.jp/~shino/toolbox/english.htm
% 
% % optN = sshist(x); hist(x,optN);
% % optW = sskernel(x); ksdensity(x,'width',optW);


% function fRates = calc_fRates(spikeTimes, nStimuli, nTransPS, nCells, transTS)
% fRates = zeros(nTransPS, nCells, nStimuli);
% edges = (0 : (nTransPS * nStimuli)) * transTS; % Calculate bin edges
% if isvector(spikeTimes) % Pad into matrix
%     spikeTimes = [spikeTimes, (ones(length(spikeTimes),1)*-1)]; % This forces histc to bin spikes for each cell even when there is at most one spike per cell
% end
% %spikeTimes(spikeTimes==0) = -1; % Don't count zeros padding the matrix
% % No longer necessary since spikes are padded with -1
% freqs = histc(spikeTimes',edges); %n(k) counts the value x(i) if edges(k) <= x(i) < edges(k+1). 
% %The last bin counts any values of x that match edges(end).
% %freqs(end-1,:)=freqs(end-1,:)+freqs(end,:); % include spikes at last timestep - not necessary since using C numbering last possible spike is edges(end)-1
% freqs(end,:) = []; % remove n^th +1 bin (spikeTimes==edges(end)) %freqs = freqs(1:end-1,:); 
% %freqs = freqs * 1/transP_Test; %-timewin makes this redundant
% for s=1:nStimuli
%     fRates(:,:,s) = freqs(1+(s-1)*nTransPS:s*nTransPS,:);
% end


function nColours = plotConnectivity(cnx, weights, MP, scaleWidth)
LWS = 4; % Line width scale
mSize = 5;
margin = 0.05;
%handle = figure();
colormap(jet);
cmp = colormap;
nColours=size(cmp,1);
if nargin < 4
    scaleWidth=false;
end

for l=0:MP.nLayers-1 % Plot excitatory neurons
    hold on;
    plot(0:MP.vExcit(l+1)-1,l,'ok','LineStyle','none','MarkerFaceColor','g','MarkerSize',mSize); % use draw?
    
    if (l>0)
        for n=1:MP.vExcit(l+1)
            nSyn = cnx{l}(n,1);
            for s=1:nSyn
                line = plot([n-1,cnx{l}(n,s+1)],[l,l-1],'Color',cmp(ceil(weights{l}(n,s)*(nColours)),:));%,'LineWidth',LWS/2); %,'LineWidth',LWS*weights{l}(n,s)); % ([Xdata],[Ydata])
                if (scaleWidth)
                    set(line,'LineWidth', LWS * weights{l}(n,s));
                end
            end
        end
    end
    hold off;
end

    %colorbar;
    set(gca,'YTick',0:MP.nLayers-1);
    xlabel('Neuron');
    ylabel('Layer');
    box on;
axis([-margin max(MP.vExcit)-1+margin -margin MP.nLayers-1+margin]);


function [xcoords,ycoords] = calcNetCoords(layDim) %, markProp)
sx=1/layDim(2); % nCols
sy=1/layDim(1); % nRows
xcoords=(sx/2)+(sx*(0:m2c(layDim(2))));
ycoords=(sy/2)+(sy*(m2c(layDim(1)):-1:0))'; % Order is reversed
if length(ycoords)==1 % Row vector layer
    ycoords=ycoords*ones(length(xcoords),1); % Column vector
elseif length(xcoords)==1
    xcoords=xcoords*ones(1,length(ycoords)); % Row vector
end
% If 1D layer: [1,vExcit(l)]
%zcoords=layer*ones(length(xcoords));
%[X,Y,Z] = meshgrid(xcoords, ycoords, zcoords);

% narray[l][n].x = (sp_x/2 + (narray[l][n].col * sp_x)) * mp->spatialScale; // col: 0,1,...,rowLen-1
% narray[l][n].y = (sp_y/2 + (((colLen - 1) - narray[l][n].row) * sp_y)) * mp->spatialScale; // row indices and y coordinates run opposite


function fh = plotNet3D(MP, weights, cnx, fh, layer, neurons)
% 3D plot: Plot all connections to/from a particular neuron
% Use C indexing i.e. pass layer and neurons starting at 0

if nargin <= 3 % Plot whole net
    fh = figure();
    wholeNet = true;
    recurse = false;
else
    figure(fh);
    wholeNet = false;
    recurse = true;
end

ticks = 0:0.25:1;
hold on
colormap(jet);
cmp = colormap;
nColours=size(cmp,1);
coords = cell(MP.nLayers,2);%3);

% Calculate coordinates and plot cell bodies  % Redo indexing!!!!!
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
        end
    end
end

if ~wholeNet
    if layer>0  % Draw all presynaptic connections
        plotPreCnx3D(MP, fh, cmp, nColours, coords, cnx, weights, layer, neurons, recurse);
    end
    if layer<MP.nLayers-1  % Draw all postsynaptic connections
        disp('Implement efferent connections!'); % To do!
    end
end

hold off
colorbar('YTick',1+ceil(nColours*ticks),'YTickLabel',ticks);
xlabel('x'); xlim([0,MP.spatialScale]);
ylabel('y'); ylim([0,MP.spatialScale]);
zlabel('Layer'); zlim([0,m2c(MP.nLayers)]);


function plotPreCnx3D(MP, fh, cmp, nColours, coords, cnx, weights, l, neurons, recurse)

figure(fh);
dim = MP.layDim;
% for lay=1:MP.nLayers
%     if MP.vSquare(lay)
%         dim{lay} = fliplr(dim{lay}); %dim{lay}(end:-1:1)
%     end
% end 
dim{~MP.vSquare} = fliplr(dim{~MP.vSquare}); %Correct indexing for 1D layers
for n=1:length(neurons)
    [c,r] = ind2sub(dim{c2m(l)},c2m(neurons(n))); % swap indices to transpose (c is row major)
    nSyn = cnx{l}(c2m(neurons(n)),1);
    preNeurons = cnx{l}(c2m(neurons(n)),1+(1:nSyn)); % List of connected presyn neurons
    for s=1:nSyn %length(preNeurons)
        preN = c2m(preNeurons(s));
        [preC,preR] = ind2sub(dim{l},preN); % swap indices to transpose (c is row major)
        line([coords{c2m(l),1}(c),coords{c2m(l-1),1}(preC)],...
            [coords{c2m(l),2}(r),coords{c2m(l-1),2}(preR)],[l,l-1],...
            'Color',cmp(ceil(weights{l}(c2m(neurons(n)),s)*(nColours)),:)); % Change s
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

% function orbit(deg)
% [az el] = view;
% rotvec = 0:deg/10:deg;
% for i = 1:length(rotvec)
%     view([az+rotvec(i) el])
%     drawnow
% end
