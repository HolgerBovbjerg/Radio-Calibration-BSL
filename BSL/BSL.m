%% MCMC-BSL ALGORITHM:
clear
load('Prior_min_max_values.mat')
load('Theta_true_values.mat')

N = 300; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
%% Find first proposed Theta and covariance for proposal distribution
[covariance, theta_curr] = find_cov_prior(prior);
theta_curr(3) = 1.76e8;
theta_start = theta_curr;
% As it has been found that the above estimated covariance is too large a
% smaller one is loaded
load('covariance_small_prior.mat')
covariance = covariance/1e4;
% A good choice covariance is very important to get good MCMC runs
%% "Observed data for testing"
load('S_obs_sim')
%%
k = 2500;    % Number of MCMC steps
L = 300;     % Numberof statistics vectors used per likelihood.

accept = 0;
s_sim = zeros(L,9);
thetas = zeros(4,k);
thetas(:,1) = theta_curr';

%% Calculate first likelihood (first point in markov chain)
for i = 1:L
    [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_curr);
    s_sim(i,:) = create_statistics(Pv, t);
end

loglikelihood = synth_loglikelihood(s_obs,s_sim);

%% Run MCMC
for j = 2:k
    cla reset
    % Draw a new proposed theta from multivariate normal distribution
    theta_prop = mvnrnd(theta_curr,covariance);
    
    % Check to see if the proposed theta is within the prior range
    % if not keep drawing new values  
    i = 0;
    while(check_params(theta_prop,prior)==2) 
        theta_prop = mvnrnd(theta_curr,covariance);
    end
    % Create statistics based on the proposed theta 
    parfor i = 1:L
        [Pv, t] = sim_turin_matrix_gpu(N, B, Ns, theta_prop);
        s_sim(i,:) = create_statistics(Pv, t);
    end
    % Calculate new log-likelihood from statistics
    loglikelihoodnew = (synth_loglikelihood(s_obs,s_sim));
    
    % Compare new and old likelihood and accept based on Hasting-ratio
    if exp(loglikelihoodnew-loglikelihood) > rand
        loglikelihood = loglikelihoodnew;
        theta_curr = theta_prop;
        accept = accept+1;
    end
    % Keep the most recent theta
    thetas(:,j) = theta_curr';
    real_time_plots(theta_true,thetas,j-1,accept,k,prior,theta_prop,loglikelihoodnew,loglikelihood);
end
