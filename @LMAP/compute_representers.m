function obj = compute_representers(obj)

% Various init variables
nobs = obj.physics_class.nobservations;
nsensors = obj.physics_class.nsensors;
times = obj.RTMflow_class.times;
notimes = length(times);

parfor k = 1:nobs * nsensors
    
    R_i = compute_representer_i(k);

end


end


function R_i = compute_representer_i(i)

end