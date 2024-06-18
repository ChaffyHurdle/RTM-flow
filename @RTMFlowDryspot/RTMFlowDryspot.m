classdef RTMFlowDryspot < RTMFlow


    properties

        %% Void tracking featrues
        num_voids = 0;
        void_volume;
        is_volume_void;
        is_dryspots; 


    end



    methods
        
        function obj = RTMFlowDryspot(Delaunay_mesh_class,physics_class,pressure_class)

            obj@RTMFlow(Delaunay_mesh_class,physics_class,pressure_class);
    
            %% Setting up void tracking properties
            obj.void_volume = zeros(obj.Delaunay_mesh_class.num_nodes,2);
            obj.is_volume_void = ones(obj.Delaunay_mesh_class.num_nodes,1);

        end


    end

end