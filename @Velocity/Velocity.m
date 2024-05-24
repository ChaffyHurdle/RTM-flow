classdef Velocity

    properties
    
        volume_class;
        darcy_class;

        %% piecewise constant velocity field
        time;
        velocity;

    end % end properties

    methods
        %% Methods of the Volume class

        function obj = Velocity(volume_class,darcy_class)

            obj.volume_class = volume_class;
            obj.darcy_class = darcy_class;
            
            num_cells = size(volume_class.element_connectivity,1);

            obj.time = 0;
            obj.velocity = zeros(num_cells,2);
        end    

    end % end methods

end