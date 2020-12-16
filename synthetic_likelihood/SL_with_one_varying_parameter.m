load('Prior_min_max_values.mat')
load('Theta_true_values.mat')
load('S_obs_sim.mat')


N = 300; % Number of Turin simulations
Ns = 801; % Number of sample points per Turin simulation
B = 4e9; % Bandwidth of signal: 4 GHz
L = 500; % Number of summaries per likelihood

%% To generate simulated observed data uncomment the follwing:
% [Pv, t] = sim_turin_matrix_gpu(1000, B, Ns, theta_true);
% s_obs = create_statistics(Pv, t);


%% Test likelihood
steps = 20;
Ts = linspace(pow2db(prior(1,1)),pow2db(prior(1,2)),steps);
G0s = linspace(prior(2,1),prior(2,2),steps);
lambdas = linspace(prior(3,1),prior(3,2),steps);
sigma_Ns = linspace(prior(4,1),prior(4,2),steps);
loglikelihood = zeros(steps,1);
s_sim = zeros(L,9);
tstart = tic;
for j = 1:steps
    %     Uncomment parameter to test likelihood over
    %     T = Ts(j);
    %     G0 = db2pow(G0s(j));
    %     lambda = lambdas(j);
    sigma_N = sigma_Ns(j);
    for i = 1:L
        %         Uncomment parameter to test likelihood over
        %         Vary T
        %         theta_prop = [T theta_true(2) theta_true(3) theta_true(4)];
        %         Vary G0
        %         theta_prop = [theta_true(1) G0 theta_true(3) theta_true(4)];
        %         Vary lambda
        %         theta_prop = [theta_true(1) theta_true(2) lambda theta_true(4)];
        %         Vary sigma_N
        theta = [theta_true(1) theta_true(2) theta_true(3) sigma_N];
        %         All true parameter
        %         theta_prop = [theta_true(1) theta_true(2) theta_true(3) theta_true(4)];
        [Pv, t] = sim_turin_matrix(N, Bw, Ns, theta);
        % For GPU acceleration
        % [Pv, t] = sim_turin_matrix_gpu(N, Bw, Ns, theta);
    end
    loglikelihood(j) = synth_loglikelihood(s_obs,s_sim);
end

%% plot
figure
% plot(Ts,loglikelihood)
% plot(G0s,loglikelihood)
% plot(lambdas,loglikelihood)
plot(sigma_Ns,loglikelihood)
xlabel('Parameter value')
ylabel('loglikelihood: p(s_y)|\theta)')