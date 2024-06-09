classdef RTMFlowDryspot


    properties

        %% Void tracking featrues
        num_voids = 0;
        void_volume;
        is_volume_void;
        is_dryspots; 


    end



    methods
        
        function obj = RTMFlowDryspot()
    
            %% Setting up void tracking properties
            obj.void_volume = zeros(mesh_class.num_nodes,2);
            obj.is_volume_void = ones(mesh_class.num_nodes,1);

        end


    end

end


%% Compute any voids/vacuum
    obj = obj.find_dryspots();
    obj = obj.void_partition();
    obj = obj.apply_ideal_gas_law();

