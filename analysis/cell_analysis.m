function cell_analysis(MP, layer, neuron, spikes, phase, loop, DT, TotalMS)

%FTHRSH_E = 0.032; % Print out from program

% global nRecordsPL
% global loops
% global nLayers
% 
% defaults;
% parameters;

% global TotalMS;
% global DT;

% if ((output_rec < 1) || (output_rec > MP.nRecordsPL))
%     disp(['Please enter an ouput cell record number from 1 to ', ...
%         int2str(MP.nRecordsPL)]);
%     exit;
% end
% 
% if ((loop < 1) || (loop > MP.loops))
%     disp(['Please enter a training loop number from 1 to ', ...
%         int2str(MP.loops)]);
%     exit;
% end

%rows = 1+(2*nRecordsPL);
%rows = 2+MP.nRecordsPL;
time_bins = 1:TotalMS; % 0?
xlimits = [1,TotalMS];


Rec = load_records(MP, layer, neuron, loop);
figure();
%subplot(4,1,1);
plot(time_bins,Rec.V(:,loop+1));
hold on

%spikes = ESpikes{layer+1}(neuron+1,:);

% if (strcmpi(phase, 'Testing'))
%     spikes = ESpikes{layer+1}(neuron+1,:);
% elseif (strcmpi(phase, 'Training'))
%     spikes = dlmread(strcat('E',int2str(loop),'L',int2str(layer),'ExcitSpikes.dat'));
%     spikes = spikes(neuron+1,2:end);
% end
%[neurons, maxspikes] = size(spikes);

s = 1;
while ((s <= length(spikes)) && (spikes(s) ~= 0))
    plot([(spikes(s)*DT*1000),(spikes(s)*DT*1000)],[MP.ThreshE-0.005,0],'-b');
    s = s + 1;
end
% http://blogs.mathworks.com/pick/2006/05/26/plotting-multiple-y-scales/

% hold all

xlabel('Times (ms)');
xlim(xlimits);
ylabel('Cell Potential');
ylim([MP.VrevI-0.005,0.005]);

hold off




% figure();
% % Plot the cell potential for a given *output layer* record and loop
% subplot(rows,1,1);
% plot(time_bins,RCELLV{nLayers,output_rec}(loop,:), 'b');
% hold on;
% % Plot cell potential threshold
% plot([0,TotalMS],[FTHRSH_E, FTHRSH_E], 'g:');
% xlabel('Times (ms)');
% xlim(xlimits);
% ylabel('Cell Potential');
% %ylim([0,1]);
% 
% % Plot learning variable 'D' for the output layer record
% subplot(rows,1,2);
% plot(time_bins,RD{nLayers,output_rec}(loop,:), 'k-');
% %raster(output_spikes,T,DT); % Training spikes are not recorded
% xlabel('Times (ms)');
% xlim(xlimits);
% ylabel('D');
% ylim([0,1]);
% 
% % For each (connected) pre-synaptic record, plot the weight and 'C'
% for r=1:nRecordsPL
% % Plot only if there is a presynaptic connection?
%     subplot(rows,1,2+r);
%     % Plot Weights
%     %NWLAYERS = nLayers-1;
% %     [row,col] = find(FLAGS);
% %     tmp = [row';col'];
% %     tmp = sortrows(tmp',1)';
% %     row = tmp(1,:);
% %     col = tmp(2,:);
% %     lay_ind = row(((NWLAYERS-1)*nRecordsPL)+r); % Previous layer's records
% %     rec_ind = col(((NWLAYERS-1)*nRecordsPL)+r);
%     lay_ind = nLayers-1; % Penultimate layer
%     neu_ind = REC_IND(lay_ind,r); % Presynaptic neuron
%     postsyn_n = REC_IND(lay_ind+1,output_rec); % Postsynaptic neuron
%     syn_ind = find(PRESYNCNX_EE(:,postsyn_n)==neu_ind); % Synaptic index
%     plot(time_bins,RWEIGHTS{lay_ind,output_rec}(syn_ind,:,loop), 'b');
%     
%     hold on;
%     % Plot C
%     plot(time_bins,RC{lay_ind,r}(loop,:), 'k-');
% 
%     xlabel({['L', int2str(lay_ind),'N',int2str(neu_ind),' -->|',int2str(syn_ind),'| --> L',...
%         int2str(lay_ind+1),'N',int2str(postsyn_n),...
%         ' | Mean=',num2str(RWSTATS{lay_ind,output_rec}(syn_ind,1,loop)),...
%         ', SD=',num2str(RWSTATS{lay_ind,output_rec}(syn_ind,2,loop))];'Times (ms)'});
%     xlim(xlimits);
%     ylabel('Presynaptic weight & C');
%     ylim([0,1]);
%     
%     %subplot(1+(2*r)); % Conductance plots
% end