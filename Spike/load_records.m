function record = load_records(layer, neuron, loops)

filename = strcat(['rL',int2str(layer),'N',int2str(neuron),'cellV.dat']);
record.cellV = load(filename);
% fid = fopen(filename, 'rb');
% tmp = fread(fid, 'float64');
% record.cellV = reshape(tmp, TotalMS, loops)';
% fclose(fid);

% if (layer<nLayers-1)
if (layer>0)
    % Read in C for presynaptic cells
    filename = strcat(['rL',int2str(layer),'N',int2str(neuron),'C.dat']);
    tmp = load(filename);
    % tmp = reshape(tmp, TotalMS, [], loops);
    block = size(tmp,1)/loops;
    for l=1:loops
        record.C(:,:,l) = tmp(1+((l-1)*block):l*block,:);
    end
    
    % Read in g for presynaptic cells
    filename = strcat(['rL',int2str(layer),'N',int2str(neuron),'g.dat']);
    tmp = load(filename);
    block = size(tmp,1)/loops;
    for l=1:loops
        record.g(:,:,l) = tmp(1+((l-1)*block):l*block,:);
    end

    % Read in synaptic weights
    filename = strcat(['rL',int2str(layer),'N',int2str(neuron),'dg.dat']);
    tmp = load(filename);
    block = size(tmp,1)/loops;
    for l=1:loops
        record.dg(:,:,l) = tmp(1+((l-1)*block):l*block,:);
    end

    % Read in D for postsynaptic cells
    filename = strcat(['rL',int2str(layer),'N',int2str(neuron),'D.dat']);
    record.D = load(filename);
end