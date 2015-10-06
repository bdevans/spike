function fh = cellAnalysis(MP, layer, neuron, spikes, prefix, loop, TotalMS, sInj, nInj)

%spikes = spikes * DT * 1000; %Outside

if nargin > 7
    plotStim = true;
else
    plotStim = false;
end

if ((loop < 0) || (loop > MP.loops-1))
    disp(['Please enter a training loop number from 0 to ', ...
        int2str(MP.loops-1)]);
    exit;
end

if strcmp(prefix(1:2),'RE') % Learning epoch
    learning = true;
else
    learning = false;
end

superimpose = true; % http://blogs.mathworks.com/pick/2006/05/26/plotting-multiple-y-scales/
ws = 2.5; % LineWidth scaling for D and \Sigma g
time_bins = 0:TotalMS; % 0?
xlimits = [0,TotalMS];
spikes(spikes<0)=[]; % Remove padding

Rec = loadRecords(MP, prefix, layer, loop); %prefix was target

ISI(spikes);
%isifh = ISI(spikes); saveas(isifh,[prefix,'ISIhist'],'png'); close(isifh);


% Plot instantaneous firing rate
% iff=zeros(length(xlimits));
% iff(1+(1:length(diff(spikes)))=diff(spikes);
%iff=[0,1000./diff(spikes)];
%if length(iff) < length(xlimits); iff(length(xlimits))=0; end % Pad with 0's
%iff=[iff,zeros(1,length(xlimits)-length(iff))];
iff = diff(spikes);
if ~isempty(iff)
    fh = figure(); stairs(spikes,[0,1000./iff]); %stairs(iff);
    hold on;
    plot(spikes,zeros(1,length(spikes)),'r*'); % Plot spikes too
    hold off;
    xlim(xlimits);          xlabel('Times (ms)');
    ylabel('Instantaneous firing frequency (Hz)');
    title(['L',int2str(layer),'N',int2str(neuron),' Spikes and IFF']);
    saveas(fh,[prefix,'IFF'],'png');    close(fh);
end
% Filter the IFF (low pass) to smooth artifacts from spikes falling either
% side of bin edges
% taus=[0.1 1.0 10 100];
% clf; hold on
% clr='grbk';
% for k=1:4
% [b, a]=butterlow1(t_samp / taus(k));
% ffr = filtfilt(b, a, ifr);
% plot(tt, ffr, clr(k));
% end
% legend('\tau = 0.1 s', '\tau = 1.0 s', '\tau = 10 s', '\tau = 100 s');

fh = figure();
tString = ['L',int2str(layer),'N',int2str(neuron)];
p=1;
if (layer > 0) % C numbering
    if learning % Learning epoch (include C, D & Dg)
        nPlots = 4;    
    else % Just plot Cell membrane potential and EfE conductances
        nPlots = 2;
    end
else % Just plot Cell membrane potential
    nPlots = 1;
end
if plotStim % Add injected current plot
    nPlots = nPlots+1;
end
if isfield(Rec,'LsigGI')
    plotSigGI = true;
    nPlots = nPlots+1;
else
    plotSigGI = false;
end
if (MP.adaptation) % Plot cell firing-rate adaptation
    nPlots = nPlots+1;
end
if MP.pCnxElE(layer+1) > eps % Plot conductances
    nPlots = nPlots+1;
    if learning && MP.trainElE % Also plot plasticity variables and weights
        nPlots = nPlots+2;
    end
end

if plotStim % Plot directly injected current
    subplot(nPlots,1,p);
    hold all
    stairs(time_bins,sInj'*10^9); %,'r--' Plot current for all groups/stimuli
    
    nStim = size(sInj,1);
    if layer==0 % Plot input neuron's injected current
        stairs(time_bins,nInj'*10^9,'r-'); % Plot current for this neuron
%         hold on;
%         stairs(time_bins,sInj'*10^9,'r--'); % Plot current for group/stimulus
%         hold off;
%     else
%         % Plot each stimulus/group
%         hold all;
%         for s=1:MP.nStimuli % Depends on training/testing
%             stairs(time_bins,sInj(s,:)*10^9);
%         end
%         hold off;
        labels = cell(1,nStim+1);
        labels{nStim+1} = tString; %['L',int2str(layer),'N',int2str(neuron)];
    else
        labels = cell(1,nStim);
    end
    hold off
    
    for s=1:nStim
        labels{s} = ['Stimulus ',int2str(s)];
    end
    legend(labels);
    xlim(xlimits);
    ylim([-0.05*MP.current*10^9,MP.current*10^9*1.05]);
    ylabel('I_{ext} (nA)');
    box on;
    if p==nPlots
        xlabel('Times (ms)');
    else
        set(gca,'XTickLabel',[]);
        p=p+1;
    end
end

% Plot V and spikes
subplot(nPlots,1,p);
plot(time_bins,Rec.V(:,loop+1)*10^3);
hold on
spikes(spikes(:)<0) = []; % Remove padding
maxspikes = length(spikes);
ypos = zeros(2,maxspikes);
ypos(1,:) = MP.ThreshE*10^3 - 5;
plot([spikes;spikes],ypos,'-b');
xlim(xlimits);
if p==nPlots %~MP.adaptation && layer==0
    xlabel('Times (ms)');
else
    set(gca,'XTickLabel',[]); 
    p=p+1;
end
ylabel('Cell Potential (mV)');
ylim([min([MP.VrevI,MP.VrestE,MP.VhyperE])*10^3-5,5]);
hold off

if (MP.adaptation) % Plot adaptation
    subplot(nPlots,1,p);
    plot(time_bins,Rec.cCa(:,loop+1)*10^6,'-g');
    xlim(xlimits);
    if layer==0 && ~MP.trainElE
        xlabel('Times (ms)');
    else
        set(gca,'XTickLabel',[]); 
    end
    ylabel('[Ca^{2+}] (\muM)');
    p=p+1;
end

if plotSigGI || layer > 0 || (learning && MP.trainElE) % Plot \Sigma g
    subplot(nPlots,1,p)
    labels = {}; %cell(3,1); 
    c=0;
    if plotSigGI % Plot sum of inhibitory conductances
        plot(time_bins,Rec.LsigGI(:,loop+1)*10^6,'-r'); % Plot in micro S
        %ylabel('$$\Sigma g_I $$ ($$\mu$$S)','interpreter','latex');
        c=c+1;
        labels{c} = '\Sigma g_I';
    end
    if layer > 0
        if (c == 1) hold all; end
        plot(time_bins,sum(Rec.Fg(:,:,loop+1))*10^6,'-b'); 
        c=c+1;
        labels{c} = '\Sigma g_{EfE}';
    end
    if MP.pCnxElE(layer+1) > eps
        if (c == 1) hold all; end
        plot(time_bins,sum(Rec.Lg(:,:,loop+1))*10^6,'-g');
        c=c+1;
        labels{c} = '\Sigma g_{ElE}';
    end
    xlim(xlimits);
    ylabel('Sum of conductances ($\mu$S)','interpreter','latex');
    if (c > 1) 
        hold off; 
        legend(labels);
    end
    if p==nPlots
        xlabel('Times (ms)');
    else
        set(gca,'XTickLabel',[]);
        p=p+1;
    end
end

if MP.pCnxElE(layer+1) > eps % Plot lateral (ElE) synaptic variables
    % Plot pre-synaptic conductances
    subplot(nPlots,1,p)
    
    if superimpose
        [AX,H1,H2] = plotyy(time_bins,Rec.Lg(:,:,loop+1)*10^9,time_bins,sum(Rec.Lg(:,:,loop+1))*10^6);
        set(AX(1),'xlim',xlimits);  set(AX(2),'xlim',xlimits);
        set(get(AX(1),'Ylabel'),'String','$g_{ElE}$ (nS)','interpreter','latex');
        set(get(AX(2),'Ylabel'),'String','$\Sigma g_{ElE} $ ($\mu$S)','interpreter','latex');
        set(H2,'LineStyle','-','LineWidth',ws*get(H1(1),'LineWidth')); %'Color',[1,1,1],
    else
        plot(time_bins,Rec.Lg(:,:,loop+1)*10^9);
        ylabel('$g_{lat}$ (nS)','interpreter','latex');
        xlim(xlimits);          %set(gca,'XTickLabel',[]);%xlabel('Times (ms)');
    end
    if p==nPlots
        xlabel('Times (ms)');
    else
        if superimpose
            set(AX(1),'XTickLabel',[]);
            set(AX(2),'XTickLabel',[]);
        else
            set(gca,'XTickLabel',[]);
        end
        p=p+1;
    end
    
    if learning && MP.trainElE
        % Plot pre-synaptic C & post-synaptic D
        subplot(nPlots,1,p)
        plot(time_bins,Rec.LC(:,:,loop+1));
        hold all
        plot(time_bins,Rec.D(:,loop+1),'--k','LineWidth',ws*get(gca,'LineWidth'));
        hold off
        xlim(xlimits);          set(gca,'XTickLabel',[]);%xlabel('Times (ms)');
        ylim([-0.05,1.05]);     ylabel('Lateral Plasticity');
        p=p+1;

        % Plot pre-synaptic weights
        subplot(nPlots,1,p)
        plot(time_bins,Rec.Ldg(:,:,loop+1));
        xlim(xlimits);          %set(gca,'XTickLabel',[]);%xlabel('Times (ms)');
        ylim([-0.05,1.05]);     ylabel('$\Delta g_{ElE}$','interpreter','latex');
        if superimpose
            hold on;
            plot(time_bins,mean(Rec.Ldg(:,:,loop+1)),'--k','LineWidth',ws*get(gca,'LineWidth'));
            hold off;
        end
        
        if p==nPlots
            xlabel('Times (ms)');
        else
            set(gca,'XTickLabel',[]);
            p=p+1;
        end
    end
end

if (layer > 0) % C numbering
    % Plot pre-synaptic EfE conductances
    subplot(nPlots,1,p)
    if superimpose
        [AX,H1,H2] = plotyy(time_bins,Rec.Fg(:,:,loop+1)*10^9,time_bins,sum(Rec.Fg(:,:,loop+1))*10^6);
        set(AX(1),'xlim',xlimits);  set(AX(2),'xlim',xlimits);
        set(get(AX(1),'Ylabel'),'String','$g_{EfE}$ (nS)','interpreter','latex');
        set(get(AX(2),'Ylabel'),'String','$\Sigma g_{EfE} $ ($\mu$S)','interpreter','latex');
        set(H2,'LineStyle','-','LineWidth',ws*get(H1(1),'LineWidth')); %'Color',[1,1,1],
    else
        plot(time_bins,Rec.Fg(:,:,loop+1)*10^9);
        ylabel('$g_{ff}$ (nS)','interpreter','latex'); 
        xlim(xlimits);          %set(gca,'XTickLabel',[]);%xlabel('Times (ms)');
    end
    if p==nPlots
        xlabel('Times (ms)');
    else
        if superimpose
            set(AX(1),'XTickLabel',[]);
            set(AX(2),'XTickLabel',[]);
        else
            set(gca,'XTickLabel',[]);
        end
        p=p+1;
    end
    
    if learning % Learning epoch
        % Plot pre-synaptic C & post-synaptic D
        subplot(nPlots,1,p)
        plot(time_bins,Rec.FC(:,:,loop+1));
        hold all
        plot(time_bins,Rec.D(:,loop+1),'--k','LineWidth',ws*get(gca,'LineWidth'));
        hold off
        xlim(xlimits);          set(gca,'XTickLabel',[]);%xlabel('Times (ms)');
        ylim([-0.05,1.05]);     ylabel('FF Plasticity');
        p=p+1;

        % Plot pre-synaptic weights 
        subplot(nPlots,1,p)
        plot(time_bins,Rec.Fdg(:,:,loop+1));
        xlim(xlimits);          xlabel('Times (ms)');
        ylim([-0.05,1.05]);     ylabel('$\Delta g_{EfE}$','interpreter','latex');
        if superimpose
            hold on;
            plot(time_bins,mean(Rec.Fdg(:,:,loop+1)),'--k','LineWidth',ws*get(gca,'LineWidth'));
            hold off;
        end
    end
    
    if p==nPlots
        xlabel('Time (ms)');
    else
        set(gca,'XTickLabel',[]);
    end
    mtit([tString,' with L',int2str(layer-1),' (afferent) variables']);
else
    mtit(prefix);%(['L',int2str(layer),'N',int2str(neuron)]);
end



function record = loadRecords(MP, prefix, layer, loop)
filename = strcat([prefix,'V.dat']);
record.V(:,loop+1) = load(filename);
% fid = fopen(filename, 'rb'); tmp = fread(fid, 'float64'); fclose(fid);
% record.cellV = reshape(tmp, TotalMS, loops)';

if strcmp(prefix(1:2),'RE') % Learning epoch
    learning = true;
else
    learning = false;
end

if MP.adaptation
    filename = strcat([prefix,'cCa.dat']);
    record.cCa(:,loop+1) = load(filename);
end

if learning
    % Read in D for postsynaptic cells
    filename = strcat([prefix,'D.dat']);
    record.D(:,loop+1) = load(filename);
end

if exist([prefix,'LAffsigGI.dat'],'file') % Load sum of inhibitory conductances
    record.LsigGI(:,loop+1) = load([prefix,'LAffsigGI.dat']);
end

if (layer>0) % Synaptic variable records associated with post-syn neuron
    % Read in g for presynaptic cells
    filename = strcat([prefix,'Affg.dat']);
    if ~exist(filename,'file'); filename=strcat([prefix,'FAffg.dat']); end
    record.Fg(:,:,loop+1) = load(filename);
    
    if strcmp(prefix(1:2),'RE') % Learning epoch
        % Read in C for presynaptic cells 
        filename = strcat([prefix,'AffC.dat']);
        if ~exist(filename,'file'); filename=strcat([prefix,'FAffC.dat']); end
        record.FC(:,:,loop+1) = load(filename); % was record.C
    %     tmp = load(filename);%     % tmp = reshape(tmp, TotalMS, [], loops);
    %     block = size(tmp,1)/loops;
    %     for l=1:loops; record.C(:,:,l) = tmp(1+((l-1)*block):l*block,:); end

        % Read in synaptic weights
        filename = strcat([prefix,'Affdg.dat']);
        if ~exist(filename,'file'); filename=strcat([prefix,'FAffdg.dat']); end
        record.Fdg(:,:,loop+1) = load(filename);
    end
end

if MP.pCnxElE(layer+1) > eps % Change to nLAffs > 0
    % Read in g for presynaptic cells
    filename = strcat([prefix,'LAffg.dat']);
    record.Lg(:,:,loop+1) = load(filename);
    
    if learning && MP.trainElE
        % Read in C for presynaptic cells 
        filename = strcat([prefix,'LAffC.dat']);
        record.LC(:,:,loop+1) = load(filename);

        % Read in synaptic weights
        filename = strcat([prefix,'LAffdg.dat']);
        record.Ldg(:,:,loop+1) = load(filename);
    end
end
