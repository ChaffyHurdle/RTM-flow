function [allSamples,timer,iterations] = run_MCMC(u0,C,params,experiment,plotting)

d = reshape(experiment.d,[],1);
Sigma = reshape(experiment.Sigma,[],1);
Sigma_minus_half = diag(1./sqrt(Sigma));
Sigma = diag(Sigma);

% Parameters
numChains = 10;             % Number of chains to run in parallel
numSamples = 100000;         % Number of MCMC samples per chain
initialGuesses = mvnrnd(u0,C,numChains); % Random initial states (d: dimension)

delete(gcp('nocreate'))
fprintf('Number of slots available: %d\n', numChains);
parpool('Threads', numChains);

% Preallocate storage for results
samples = cell(numChains, 1);

tic;

% Parallel execution of chains
parfor c = 1:numChains
    rng(c); % Set random seed for reproducibility in each chain
    chain = zeros(params.Nx,numSamples); % Preallocate for the chain
    currentSample = initialGuesses(c, :); % Initialize state
    [currentp,~] = forward_map(currentSample,params);
    currentJ = 0.5*norm(Sigma_minus_half*(d - reshape(currentp,[],1)))^2;
    xi_vec = mvnrnd(u0,C,numSamples);

    beta = 0.1;
    target_accept_rate = 0.3;   % Target acceptance rate
    adapt_interval = 100;        % Interval for adapting beta
    adapt_factor = 1.1;         % Factor for increasing/decreasing beta
    accept_count = 0;           % Count accepted proposals
    
    for n = 1:numSamples
        % Propose using pCN formula
        proposal = sqrt(1-beta^2)*currentSample + (1-sqrt(1-beta^2))*u0 + beta*xi_vec(n,:);
        [proposalp,~] = forward_map(proposal,params);
        proposalJ = 0.5*norm(Sigma_minus_half*(d - reshape(proposalp,[],1)))^2;

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
    
    samples{c} = chain(:,numSamples/10:end); % Save chain results
end

% Combine results from all chains
allSamples = horzcat(samples{:});

timer = toc

figure(4)
U_mean=mean(allSamples,2);
X=[params.x_locations,fliplr(params.x_locations)];
Y=[(U_mean-0.674*sqrt(var(allSamples,0,2)))',fliplr((U_mean+0.674*sqrt(var(allSamples,0,2)))')];
fill(X,Y, 'g', 'EdgeColor','none', 'FaceAlpha',0.25);
hold on
plot(params.x_locations,U_mean,'k--')
plot(experiment.params_fwd.x_locations,experiment.u_true,'r')
plot(params.x_locations,U_mean+1.96*sqrt(var(allSamples,0,2)),'k')
plot(params.x_locations,U_mean-1.96*sqrt(var(allSamples,0,2)),'k')
xline(experiment.ups_true(end),"r--")
xlim([0,1])
ylim([-2,2])
hold off
drawnow

iterations = numSamples*numChains;

end