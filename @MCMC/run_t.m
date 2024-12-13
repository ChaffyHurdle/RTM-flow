function obj = run_t(obj,t)

mesh_class = obj.mesh_class;
physics_class = obj.physics_class;
inverse_class = obj.inverse_class;

data = reshape(inverse_class.data(:,1:t),[],1);
Sigma = reshape(inverse_class.Sigma(:,1:t),[],1);
Sigma_minus_half = diag(1./sqrt(Sigma));
Sigma = diag(Sigma);

% Parameters
u0 = inverse_class.u0;
C = inverse_class.C0_inv;
nchains = obj.nchains;
nsamples = obj.nsamples;
initial_samples = mvnrnd(u0,C,nchains);

delete(gcp('nocreate'))
fprintf('Number of slots available: %d\n', nchains);
parpool('Threads', nchains);

% Preallocate storage for results
samples = cell(nchains, 1);

tic;

% Parallel execution of chains
parfor c = 1:nchains
    rng(c); % Set random seed for reproducibility in each chain
    chain = zeros(mesh_class.num_elements,nsamples); % Preallocate for the chain
    currentSample = initial_samples(c, :); % Initialize state

    physics = physics_class;
    physics.permeability = exp(currentSample)';
    pressure = Pressure(mesh_class,physics);
    RTMflow = RTMFlow(mesh_class,physics,pressure);
    RTMflow = RTMflow.run(physics.observation_times(t));

    currentJ = 0.5*norm(Sigma_minus_half*(data - reshape(RTMflow.pressure_data(:,1:t),[],1)))^2;
    xi_vec = mvnrnd(u0,C,nsamples);

    beta = 0.1;
    target_accept_rate = 0.3;   % Target acceptance rate
    adapt_interval = 100;        % Interval for adapting beta
    adapt_factor = 1.1;         % Factor for increasing/decreasing beta
    accept_count = 0;           % Count accepted proposals
    
    for n = 1:nsamples
        disp(n)
        % Propose using pCN formula
        proposal = sqrt(1-beta^2)*currentSample + (1-sqrt(1-beta^2))*u0 + beta*xi_vec(n,:);

        physics = physics_class;
        physics.permeability = exp(proposal)';
        pressure = Pressure(mesh_class,physics);
        RTMflow = RTMFlow(mesh_class,physics,pressure);
        RTMflow = RTMflow.run(physics.observation_times(t));
        proposalJ = 0.5*norm(Sigma_minus_half*(data - reshape(RTMflow.pressure_data(:,1:t),[],1)))^2;

        % Compute acceptance probability (posterior ratio)
        alpha = min(1, exp(currentJ - proposalJ));
        
        % Accept/reject step
        if rand < alpha
            currentSample = proposal; % Accept
            currentJ = proposalJ;
            accept_count = accept_count + 1;
        end
        
        chain(:,n) = currentSample; % Store sample

        if mod(n, adapt_interval) == 0
            accept_rate = accept_count / adapt_interval; % Compute acceptance rate
            if accept_rate > target_accept_rate
                beta = beta * adapt_factor; % Increase beta
            elseif accept_rate < target_accept_rate
                beta = beta / adapt_factor; % Decrease beta
            end
                accept_count = 0; % Reset accept count
        end
    end
    
    samples{c} = chain(:,nsamples/10:end); % Save chain results
end

% Combine results from all chains
allSamples = horzcat(samples{:});

timer = toc

U_mean=mean(allSamples,2);
U_var = var(allSamples,0,2);

figure(1)
subplot(1,3,1)
pdeplot(inverse_class.fwd_mesh.nodes',inverse_class.fwd_mesh.elements',XYData=inverse_class.u_true, ...
    XYStyle='interp',ColorMap="jet",Mesh="off")
clim([-1.5,1.5])
subplot(1,3,2)
pdeplot(mesh_class.nodes',mesh_class.elements',XYData=U_mean, ...
    XYStyle='interp',ColorMap="jet",Mesh="off")
clim([-1.5,1.5])
subplot(1,3,3)
pdeplot(mesh_class.nodes',mesh_class.elements',XYData=U_var, ...
    XYStyle='interp',ColorMap="jet",Mesh="off")
clim([0,0.25])
drawnow

iterations = nsamples*nchains;

end