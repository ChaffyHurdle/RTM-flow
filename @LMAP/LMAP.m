classdef LMAP

    properties
        
        % Properties associated with pressure data collection
        mesh_class;
        physics_class;
        RTMflow_true;
        inverse_class;

        u_true;
        u0;
        pressure_data;
        sigma_delta;
        delta;

    end

    methods

        function obj = LMAP(RTMflow_class, inverse_class)
            
            obj.mesh_class = RTMflow_class.Delaunay_mesh_class;
            obj.physics_class = RTMflow_class.physics_class;
            obj.RTMflow_true = RTMflow_class;
            obj.inverse_class = inverse_class;

            obj.u_true = inverse_class.u_meshcenters;
            obj.u0 = zeros(1,length(obj.u_true));
            obj.pressure_data = RTMflow_class.pressure_data;
            obj.sigma_delta = 1/100;
            obj.delta = @(x_i,x) exp(-0.5*norm(x-x_i)^2/(2*obj.sigma_delta^2))/sqrt(2*pi*obj.sigma_delta^2);

        end
    end

end