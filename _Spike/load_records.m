function record = load_records(MP, layer, neuron, loop)



% if MP.pretrain
%     prefix = strcat(['RptL',int2str(layer),'N',int2str(neuron)]);
% end
% 
% if MP.train
%     prefix = strcat(['RE',int2str(loop),'L',int2str(layer),'N',int2str(neuron)]);
% end

prefix = strcat(['Records/RE',int2str(loop),'L',int2str(layer),'N',int2str(neuron)]);

filename = strcat([prefix,'V.dat']);
record.V(:,loop+1) = load(filename);
% fid = fopen(filename, 'rb');
% tmp = fread(fid, 'float64');
% record.cellV = reshape(tmp, TotalMS, loops)';
% fclose(fid);

if MP.adaptation
    filename = strcat([prefix,'cCa.dat']);
    record.cCa(:,loop+1) = load(filename);
end

% Read in D for postsynaptic cells
filename = strcat([prefix,'D.dat']);
record.D(:,loop+1) = load(filename);


if (layer>0) % Synaptic variable records associated with post-syn neuron
    
    % Read in C for presynaptic cells
    filename = strcat([prefix,'AffC.dat']);
    record.C(:,:,loop+1) = load(filename);
%     tmp = load(filename);
%     % tmp = reshape(tmp, TotalMS, [], loops);
%     block = size(tmp,1)/loops;
%     for l=1:loops
%         record.C(:,:,l) = tmp(1+((l-1)*block):l*block,:);
%     end
    
    % Read in g for presynaptic cells
    filename = strcat([prefix,'Affg.dat']);
    record.g(:,:,loop+1) = load(filename);
    
    % Read in synaptic weights
    filename = strcat([prefix,'Affdg.dat']);
    record.dg(:,:,loop+1) = load(filename);
end