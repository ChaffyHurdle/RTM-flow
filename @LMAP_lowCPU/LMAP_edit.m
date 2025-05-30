classdef LMAP_edit

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
        u0;
        lambdas;
        grad_lambdas;
        bold_R;
        curly_R
        d;
        i_vec;
        j_vec;

        % LMAP parameters
        tol1;
        tol2;
        max_iterations;
        alpha;
        scale;

        % Saves
        u_map;
        C_map;
        u_iterations;
        J_iterations;
        scaled_data_misfit;
        execution_times;
        best_alpha;
        umap_seq;
        Cmap_seq;
        timer_seq;

    end

    methods

        function obj = LMAP_edit(inverse_class,physics_class,alpha0,scale,tolU,tolJ)
            
            % Forward simulation
            obj.inverse_class = inverse_class;
            obj.mesh_class = inverse_class.inv_mesh;
            obj.physics_class = physics_class;
            obj.physics_class.permeability = exp(inverse_class.u0)';
            obj.pressure_class = Pressure(obj.mesh_class,obj.physics_class);
            RTMflow_class = RTMFlow(obj.mesh_class,obj.physics_class,obj.pressure_class);
            RTMflow_class = RTMflow_class.run(inf);
            obj.RTMflow_class = RTMflow_class;
            
            % Dirac delta
            obj.sigma_delta_t = 0.0018/1000;
            obj.delta_t = @(t_i,t,delt_t) exp(-abs(t-t_i).^2/(2*delt_t))/sqrt(2*pi*delt_t);
           
            % LMAP
            [i_vec, j_vec] = meshgrid(1:physics_class.nsensors, 1:physics_class.nobservations);
            obj.i_vec = reshape(i_vec', [], 1);
            obj.j_vec = reshape(j_vec', [], 1);
            
            obj.u = inverse_class.u0;
            obj.u0 = inverse_class.u0;
            obj.alpha = alpha0;
            obj.scale = scale;
            obj.tol1 = tolU;
            obj.tol2 = tolJ;
            obj.max_iterations = 50;
        end
    end

end