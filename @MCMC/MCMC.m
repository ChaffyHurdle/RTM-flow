classdef MCMC

    properties
        
        % Forward simulation
        inverse_class
        mesh_class;
        physics_class;

        % MCMC parameters
        beta0;
        nsamples;
        nchains;
        target_accept_rate;
        adapt_interval;
        adapt_factor;

        % Saves
        u_mcmc;
        C_mcmc;
        u_iterations;
        J_iterations;
        execution_times;

        umcmc_seq;
        Cmcmc_seq;
        timer_seq;

    end

    methods

        function obj = MCMC(inverse_class,physics_class,nsamples,nchains)
            
            % Forward simulation
            obj.inverse_class = inverse_class;
            obj.mesh_class = inverse_class.inv_mesh;
            obj.physics_class = physics_class;
           
            % MCMC
            obj.beta0 = 0.1;
            obj.nsamples = nsamples;
            obj.nchains = nchains;
            obj.target_accept_rate = 0.3;
            obj.adapt_interval = 100;
            obj.adapt_factor = 1.1;

        end
    end

end