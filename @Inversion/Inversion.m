classdef Inversion

    properties

        % Classes
        fwd_mesh;
        inv_mesh;
        
        % Related to Matern prior used to generate K (or u = log K)
        matern_var;
        matern_length_scale;
        matern_nu;

        C0_fwd;
        cholesky_L_fwd;
        u_true;

        C0_inv;
        cholesky_L_inv;
        u0;

        pressure_noiseless;
        Sigma;
        data;
        
    end

    methods
        function obj = Inversion(fwd_mesh,inv_mesh,matern_args)
            
            % Classes
            obj.fwd_mesh = fwd_mesh;
            obj.inv_mesh = inv_mesh;
            
            % Instantiate Matern prior parameters
            obj.matern_var = matern_args(1);
            obj.matern_length_scale = matern_args(2);
            obj.matern_nu = matern_args(3);

            % Instantiate u_true, u0
            obj.u_true = zeros(1,length(fwd_mesh.centroids));
            obj.u0 = zeros(1,length(inv_mesh.centroids));

            % Generate covariance matrix from Matern parameters
            obj.C0_fwd = obj.compute_matern_covariance(fwd_mesh.centroids);
            obj.cholesky_L_fwd = chol(obj.C0_fwd,'lower');
            obj.C0_inv = obj.compute_matern_covariance(inv_mesh.centroids);
            obj.cholesky_L_inv = chol(obj.C0_inv,'lower');
        end
    end
end