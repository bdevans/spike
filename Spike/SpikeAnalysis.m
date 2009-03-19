function SpikeAnalysis(sim)

rdir = '~/Code/Spike/Results/';

path(path, '~/Code/Spike/');
eval(['cd ',strcat(rdir,sim)]); %eval(['! cd ',strcat(rdir,sim)]);

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

load L0affNeuronsEI.dat;
load L0affNeuronsElE.dat;
load L0affNeuronsIE.dat;
load L0affNeuronsII.dat;

load L1affNeuronsEI.dat;
load L1affNeuronsElE.dat;
load L1affNeuronsIE.dat;
load L1affNeuronsII.dat;
load L1affNeuronsEfE.dat;

load L1weightsEfE.dat;

L0ExcitSpikes = dlmread('L0ExcitSpikes.dat');
L1ExcitSpikes = dlmread('L1ExcitSpikes.dat');
L0InhibSpikes = dlmread('L0InhibSpikes.dat');
L1InhibSpikes = dlmread('L1InhibSpikes.dat');


% Calculate bins and centres
bins = 0.025:0.05:0.975;
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

% Plot pre-training results
if pretrain==1
    pt_output_spikes = dlmread('ptL1ExcitSpikes.dat');
    load ptL1weightsEfE.dat;
    subplot(2,3,1);
    imagesc(sort_matrix(ptL1weightsEfE'));
    xlabel('Output neuron');
    ylabel('Input neuron');

    subplot(2,3,4);
    hist(reshape(ptL1weightsEfE,1,[]),bins);
    xlabel('Synaptic weight');
    ylabel('Frequency');

    subplot(2,3,[2 6]);
    raster(pt_output_spikes,testPeriodMS,DT);
    saveas(gcf,strcat('PT',sim),'png');
end

% Plot post-training results

figure();
raster(L0ExcitSpikes,testPeriodMS,DT);
saveas(gcf,strcat('inputs',sim),'png');

figure();
subplot(2,3,1);
imagesc(sort_matrix(L1weightsEfE'));
xlabel('Output neuron');
ylabel('Input neuron');

subplot(2,3,4);
hist(reshape(L1weightsEfE,1,[]),bins);
xlabel('Synaptic weight');
ylabel('Frequency');

subplot(2,3,[2 6]);
raster(L1ExcitSpikes,testPeriodMS,DT);
saveas(gcf, sim, 'png');
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
%%
clear remain
clear token

save matlab_workspace

%%
% for f=1:NRECORDS_PL
%     cell_analysis(f,1);
% end