function obj = find_dryspots(obj)

%% unpacking variables
fFactor = obj.volume_fill_percentage;
vent_idx = obj.pressure_class.vent_idx;
nb_nodes = obj.volume_class.nb_nodes;

nnode = length(fFactor);
queue   = zeros(nnode,1);
isInQueue = false(nnode,1);
is_volume_void = zeros(nnode,1);
 
is_volume_void(fFactor == 1) = 0.5; % Filled volume

lastIdx = 0;

for i = 1 : length(vent_idx)
    vi = vent_idx(i);
    if fFactor(vi) < 1
        lastIdx = lastIdx + 1;
        queue(lastIdx) = vi;
        isInQueue(vi) = true;
    end
end


%% 

if lastIdx > 0
    startIdx = 1;
    while 1
        qs = queue(startIdx);
        
        nbqs =nb_nodes{qs};
        for i = 1 : length(nbqs)
            nbi = nbqs(i);
            if fFactor(nbi) < 1 && ~isInQueue(nbi)
                lastIdx = lastIdx+1;
                queue(lastIdx) = nbi;
                isInQueue(nbi) = true;
            end
        end
        startIdx = startIdx + 1;
        if startIdx > lastIdx
            break;
        end 
    end
    is_volume_void(queue(1:lastIdx)) = 1;
end

dryspot_flag = ~isempty(find(is_volume_void==0));

%% repacking variables
obj.is_volume_void = is_volume_void;
obj.is_dryspots = dryspot_flag;
end