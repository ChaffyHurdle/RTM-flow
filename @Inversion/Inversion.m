classdef Inversion

    properties

        % Classes
        DelaunayMesh;
        
        % Related to Matern prior used to generate K (or u = log K)
        matern_var;
        matern_length_scale;
        matern_nu;
        high_ref;
        C0;
        cholesky_L;
        u_highdim;
        u_meshcenters;
        
    end

    methods
        function obj = Inversion(DelaunayMesh,matern_args)
            
            % Classes
            obj.DelaunayMesh = DelaunayMesh;
            
            % Instantiate Matern prior parameters
            obj.matern_var = matern_args(1);
            obj.matern_length_scale = matern_args(2);
            obj.matern_nu = matern_args(3);
            obj.high_ref = 75;

            % Generate covariance matrix from Matern parameters
            obj = obj.compute_matern_covariance();
            obj.cholesky_L = chol(obj.C0,'lower');
        end
    end
end