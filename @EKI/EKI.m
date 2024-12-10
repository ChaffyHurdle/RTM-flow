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
        ensemble;
        max_iterations;
        i_vec;
        j_vec;

        u_iterations;
        J_iterations;
        scaled_data_misfit;
        execution_times;

        ueki_seq;
        Ceki_seq;
        timer_seq;

    end

    methods

        function obj = EKI(inverse_class,physics_class)
            
            % Forward simulation
            obj.inverse_class = inverse_class;
            obj.mesh_class = inverse_class.inv_mesh;
            obj.physics_class = physics_class;
           
            % EKI
            [i_vec, j_vec] = meshgrid(1:physics_class.nsensors, 1:physics_class.nobservations);
            obj.i_vec = reshape(i_vec', [], 1);
            obj.j_vec = reshape(j_vec', [], 1);
            obj.max_iterations = 50;
            obj.J = 1000;
        end
    end

end