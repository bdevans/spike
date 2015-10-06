output_rec=1;
%lay_ind=1;
loop=1;
time_bins = 1:TOT_MS;
lay_ind = NLAYERS-1; % Penultimate layer

% WSTATS=zeros(NEXCIT,2);
% WSTATS(:,1)=mean(RWEIGHTS{lay_ind,output_rec}(:,:,loop),2);
% WSTATS(:,2)=std(RWEIGHTS{lay_ind,output_rec}(:,:,loop),0,2);

postsyn_n = REC_IND(lay_ind+1,output_rec); % Postsynaptic neuron

for f=1:NEXCIT/10
    figure()
    for n=1:NEXCIT/12
        subplot(NEXCIT/12,1,n);
        %lay_ind=NLAYERS-1;
        neu_ind=n+(f-1)*(NEXCIT/12);
            %neu_ind = REC_IND(lay_ind,r); % Presynaptic neuron
        postsyn_n = REC_IND(lay_ind+1,output_rec); % Postsynaptic neuron
        syn_ind = find(PRESYNCNX_EE(:,postsyn_n)==neu_ind); % Synaptic index
        plot(time_bins,RCELLV{lay_ind+1,output_rec}(loop,:), 'b');
        hold on;
        plot(time_bins,RWEIGHTS{lay_ind,output_rec}(syn_ind,:,loop), 'k');
        xlimits = [1,TOT_MS];
        %xlabel('Times (ms)');
            xlabel({['L', int2str(lay_ind),'N',int2str(neu_ind),' -->|',int2str(syn_ind),'| --> L',...
        int2str(lay_ind+1),'N',int2str(postsyn_n),...
        ' | Mean=',num2str(RWSTATS{lay_ind,output_rec}(syn_ind,1,loop)),...
        ', SD=',num2str(RWSTATS{lay_ind,output_rec}(syn_ind,2,loop))];'Times (ms)'});
%         xlabel({['L', int2str(lay_ind),'N',int2str(neu_ind),' --> L',...
%         int2str(lay_ind+1),'N',int2str(postsyn_n),...
%         '| Mean=',num2str(RWSTATS{lay_ind,output_rec}(syn_ind,1,loop)),...
%         ', SD=',num2str(RWSTATS{lay_ind,output_rec}(syn_ind,2,loop))];'Times (ms)'});
        xlim(xlimits);
        ylabel('V and DeltaG');
        ylim([0,1]);
    end
end