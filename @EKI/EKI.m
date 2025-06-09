classdef EKI

    properties
        
        % Forward simulation
        inverse_class
        mesh_class;
        physics_class;
        
        % EKI algorithm
        u;
        u0;
        J;
        max_iterations;
        ensemble;

        u_iterations;
        J_iterations;
        scaled_data_misfit;
        execution_times;

        ueki;
        Ceki;
        timer;
        F_evals;

        ueki_seq;
        Ceki_seq;
        timer_seq;

    end

    methods

        function obj = EKI(inverse_class,physics_class,J)
            
            % Forward simulation
            obj.inverse_class = inverse_class;
            obj.mesh_class = inverse_class.inv_mesh;
            obj.physics_class = physics_class;
           
            % EKI
            obj.max_iterations = 50;
            obj.J = J;
        end
    end

end