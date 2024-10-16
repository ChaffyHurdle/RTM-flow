classdef LMAP

    properties
        
        % Forward simulation
        inverse_class
        mesh_class;
        physics_class;
        pressure_class;
        RTMflow_class;

        % Dirac delta
        sigma_delta_t;
        delta_t;
        
        % LMAP algorithm
        u;
        u0
        alpha;
        max_iterations;
        tol1;
        tol2;
        lambdas;
        grad_lambdas;
        R;
        Q;
        p_tildes;
        p_tilde0;
        i_vec;
        j_vec;
        tildePmat;
        d;
        u_map;
        C_map;
        u_iterations;
        J_iterations;
        execution_times;
        best_alpha;

    end

    methods

        function obj = LMAP(inverse_class,physics_class)
            
            % Forward simulation
            obj.inverse_class = inverse_class;
            obj.mesh_class = inverse_class.inv_mesh;
            obj.physics_class = physics_class;
            obj.physics_class.permeability = exp(inverse_class.u0)';
            obj.pressure_class = Pressure(obj.mesh_class,obj.physics_class);
            RTMflow_class = RTMFlow(obj.mesh_class,obj.physics_class,obj.pressure_class);
            RTMflow_class = RTMflow_class.run();
            obj.RTMflow_class = RTMflow_class;
            
            % Dirac delta
            obj.sigma_delta_t = 0.0018/1000;
            obj.delta_t = @(t_i,t,delt_t) exp(-abs(t-t_i)^2/(2*delt_t))/sqrt(2*pi*delt_t);
           
            % LMAP
            [i_vec, j_vec] = meshgrid(1:physics_class.nsensors, 1:physics_class.nobservations);
            obj.i_vec = reshape(i_vec', [], 1);
            obj.j_vec = reshape(j_vec', [], 1);
            obj.u = inverse_class.u0;
            obj.u0 = inverse_class.u0;
            obj.alpha = 1e5;
            obj.tol1 = 0.03;
            obj.tol2 = 0.03;
            obj.max_iterations = 50;
        end
    end

end