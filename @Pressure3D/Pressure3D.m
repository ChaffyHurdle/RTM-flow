classdef Pressure3D < Pressure

% A derived class Pressure3D from Pressure which extends the 2D FEM solver
% for the pressure field to 3D tetrahedral meshes using linear finite
% elements.

    methods


        function obj = Pressure3D(mesh_class,physics_class,inlet_func,vent_func,p_D)

            obj.mesh_class = mesh_class;
            obj.physics_class = physics_class;
            obj.inlet_func = inlet_func;
            obj.vent_func = vent_func;
            obj.p_D = p_D;

            %% set time to zero
            obj.time = 0.0;

            %% Computing information for boundary conditions
            obj = obj.compute_inlets_outlets();

            %% Allocating pressure and pressure gradient
            num_dofs = mesh_class.num_nodes;
            obj.pressure = zeros(num_dofs,1);
            obj.pressure = obj.p_D(obj);
            obj = obj.compute_shape_fun_gradients();
            obj = obj.compute_pressure_gradient();

            %% Allocating FEM system of equations
            obj.stiffness_matrix = spalloc(num_dofs,num_dofs,10*num_dofs);
            obj.load_vector = zeros(num_dofs,1);


        end
    end
end