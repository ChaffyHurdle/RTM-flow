classdef LMAP

    properties
        
        % Properties associated with pressure data collection
        inverse_class
        mesh_class;
        physics_class;
        pressure_class;
        RTMflow_class;

        sigma_delta_x;
        sigma_delta_t;
        delta_t;
        delta_x;
        f;
        
        u;
        lambdas;
        grad_lambdas;
        representers;

    end

    methods

        function obj = LMAP(inverse_class,physics_class)
            
            % Classes
            obj.inverse_class = inverse_class;
            obj.mesh_class = inverse_class.inv_mesh;
            obj.physics_class = physics_class;
            obj.physics_class.permeability = exp(inverse_class.u0);
            obj.pressure_class = Pressure(obj.mesh_class,obj.physics_class);
            obj.RTMflow_class = RTMFlow(obj.mesh_class,obj.physics_class,obj.pressure_class);
            
            % Dirac deltas
            obj.sigma_delta_t = 1e-11;
            obj.sigma_delta_x = 1/1000;
            obj.delta_t = @(t_i,t,delt_t) exp(-abs(t-t_i)^2/(2*delt_t^2))/sqrt(2*pi*delt_t^2);
            obj.delta_x = @(x_i,x) exp(-vecnorm( (x-x_i)' ).^2/(2*obj.sigma_delta_x))/sqrt(2*pi*obj.sigma_delta_x);

            % Placeholder for time-dependent boundary condition
            obj.f = @(t) sin(t/obj.RTMflow_class.T*pi)*exp(-t/obj.RTMflow_class.T);
    
        end
    end

end