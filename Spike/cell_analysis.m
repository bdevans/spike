function cell_analysis(output_rec, loop)

%FTHRSH_E = 0.032; % Print out from program

global FTHRSH_E
global PRESYNCNX_EE
global RCELLV
global RWEIGHTS
global RC
global RD
global NRECORDS_PL
global TOT_MS
global LOOPS
global NLAYERS
%global NWLAYERS
%global FLAGS
global REC_IND
global RWSTATS

if ((output_rec < 1) || (output_rec > NRECORDS_PL))
    disp(['Please enter an ouput cell record number from 1 to ', ...
        int2str(NRECORDS_PL)]);
    exit;
end

if ((loop < 1) || (loop > LOOPS))
    disp(['Please enter a training loop number from 1 to ', ...
        int2str(LOOPS)]);
    exit;
end

%rows = 1+(2*NRECORDS_PL);
rows = 2+NRECORDS_PL;
time_bins = 1:TOT_MS; % 0?
xlimits = [1,TOT_MS];




figure();
% Plot the cell potential for a given *output layer* record and loop
subplot(rows,1,1);
plot(time_bins,RCELLV{NLAYERS,output_rec}(loop,:), 'b');
hold on;
% Plot cell potential threshold
plot([0,TOT_MS],[FTHRSH_E, FTHRSH_E], 'g:');
xlabel('Times (ms)');
xlim(xlimits);
ylabel('Cell Potential');
%ylim([0,1]);

% Plot learning variable 'D' for the output layer record
subplot(rows,1,2);
plot(time_bins,RD{NLAYERS,output_rec}(loop,:), 'k-');
%raster(output_spikes,T,DT); % Training spikes are not recorded
xlabel('Times (ms)');
xlim(xlimits);
ylabel('D');
ylim([0,1]);

% For each (connected) pre-synaptic record, plot the weight and 'C'
for r=1:NRECORDS_PL
% Plot only if there is a presynaptic connection?
    subplot(rows,1,2+r);
    % Plot Weights
    %NWLAYERS = NLAYERS-1;
%     [row,col] = find(FLAGS);
%     tmp = [row';col'];
%     tmp = sortrows(tmp',1)';
%     row = tmp(1,:);
%     col = tmp(2,:);
%     lay_ind = row(((NWLAYERS-1)*NRECORDS_PL)+r); % Previous layer's records
%     rec_ind = col(((NWLAYERS-1)*NRECORDS_PL)+r);
    lay_ind = NLAYERS-1; % Penultimate layer
    neu_ind = REC_IND(lay_ind,r); % Presynaptic neuron
    postsyn_n = REC_IND(lay_ind+1,output_rec); % Postsynaptic neuron
    syn_ind = find(PRESYNCNX_EE(:,postsyn_n)==neu_ind); % Synaptic index
    plot(time_bins,RWEIGHTS{lay_ind,output_rec}(syn_ind,:,loop), 'b');
    
    hold on;
    % Plot C
    plot(time_bins,RC{lay_ind,r}(loop,:), 'k-');

    xlabel({['L', int2str(lay_ind),'N',int2str(neu_ind),' -->|',int2str(syn_ind),'| --> L',...
        int2str(lay_ind+1),'N',int2str(postsyn_n),...
        ' | Mean=',num2str(RWSTATS{lay_ind,output_rec}(syn_ind,1,loop)),...
        ', SD=',num2str(RWSTATS{lay_ind,output_rec}(syn_ind,2,loop))];'Times (ms)'});
    xlim(xlimits);
    ylabel('Presynaptic weight & C');
    ylim([0,1]);
    
    %subplot(1+(2*r)); % Conductance plots
end