classdef LMAP

    properties
        
        % Properties associated with pressure data collection
        mesh_class;
        physics_class;
        RTMflow_class;
        inverse_class;

        u_true;
        u0;
        pressure_data;
        sigma_delta_x;
        sigma_delta_t;
        delta_t;
        delta_x;
        f;

        lambdas;

    end

    methods

        function obj = LMAP(RTMflow_class, inverse_class)
            
            obj.mesh_class = RTMflow_class.Delaunay_mesh_class;
            obj.physics_class = RTMflow_class.physics_class;
            obj.RTMflow_class = RTMflow_class;
            obj.inverse_class = inverse_class;

            obj.u_true = inverse_class.u_meshcenters;
            obj.u0 = zeros(1,length(obj.u_true));
            obj.pressure_data = RTMflow_class.pressure_data;
            obj.sigma_delta_t = 1e-9;
            obj.sigma_delta_x = 1/5000;
            obj.delta_t = @(t_i,t) exp(-0.5*abs(t-t_i)^2/(2*obj.sigma_delta_t^2))/sqrt(2*pi*obj.sigma_delta_t^2);
            obj.delta_x = @(x_i,x) exp(-0.5*vecnorm( (x-x_i)' ).^2/(2*obj.sigma_delta_x^2))/sqrt(2*pi*obj.sigma_delta_x^2);
            obj.f = @(t) sin(t/obj.RTMflow_class.T*pi)*exp(-t/obj.RTMflow_class.T);

            obj.lambdas = cell(RTMflow_class.nsensors,RTMflow_class.nobservations);
        end
    end

end