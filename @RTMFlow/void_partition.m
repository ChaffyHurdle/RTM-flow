function obj = void_partition(obj)

%% skip computation if no dryspots found
if ~obj.is_dryspots
    return
end

%% unpackaging properties
num_voids = obj.num_voids;
void_volume = obj.void_volume;
fFactor = obj.volume_fill_percentage;
volume_measures = obj.volume_class.volume_measures;
Cnode = obj.volume_class.node_connectivity;
is_volume_void = obj.is_volume_void;


candidate = find(is_volume_void == 0);
biggestIdx = max(candidate);
n = length(candidate);
idx = 1;

while 1
    ci = candidate(idx);
    % new void occurs
    if void_volume(ci,1) == 0
        num_voids = num_voids + 1;
        queue = zeros(n,1);
        isInQueue = false(biggestIdx,1);
        qIdx = 1;
        volume = volume_measures(ci)*(1-fFactor(ci));
        queue(qIdx) = ci;
        isInQueue(ci) = true;
        nb_nodes = find(Cnode(ci,:));
        for j = 1 : length(nb_nodes)
            nb =nb_nodes(j);
            if is_volume_void(nb) == 0 && ~isInQueue(nb)
                qIdx = qIdx + 1;
                queue(qIdx) = nb;
                isInQueue(nb) = true;
                volume = volume + volume_measures(nb)*(1-fFactor(nb));
            end
        end
        if qIdx > 1
            startIdx = 2;
            while 1
                qs = queue(startIdx);
                nb_nodes = find(Cnode(qs,:));
                for i = 1 : length(nb_nodes)
                    nb = nb_nodes(i);
                    if is_volume_void(nb) == 0 && ~isInQueue(nb)
                        qIdx = qIdx + 1;
                        queue(qIdx) = nb;
                        isInQueue(nb) = 1;
                        volume = volume + volume_measures(nb)*(1-fFactor(nb));
                    end
                end
                startIdx = startIdx + 1;
                if startIdx > qIdx
                    break;
                end
                    
            end
        end
        void_volume(queue(1:qIdx),1) = volume;
        void_volume(queue(1:qIdx),2) = num_voids;
    end
    idx = idx + 1;
    if idx > n
        break
    end
end

%% repackaging
obj.void_volume = void_volume;
obj.num_voids = num_voids;

end