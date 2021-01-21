%% PMC-ABC ALGORITHM:
% As described in:
% Bharti, A., & Pedersen, T. (2020). Calibration of Stochastic Channel Models using Approximate Bayesian Computation.
% https://doi.org/10.1109/GCWkshps45667.2019.9024563
% (without regression step)
clear
%% --- Global turin model simulation parameters ---------------------------------------------------
N  = 300;   % Number of different turin simulations. (Corresponds to the number of columns in Pv matrix)
Ns = 801;   % Number of time entries for each turin simulation. (Corresponds to the number of rows in Pv matrix)
Bw = 4e9;   % Bandwidth (4Ghz).

%% --- Generate "observed data" -----------------------------------------------------
load("Theta_true_values.mat")

% For custom values:
% T       = 7.8e-9;
% G0      = db2pow(-83.9);
% lambda  = 10e9;
% sigma_N = 1.673e-4;
% M = 2000; % Number of summary statisctics realisations

% Theta_true_values = [T G0 lambda sigma_N];

% Load summary statistics
load('S_obs_sim_ABC.mat')
% For generating new S_obs uncomment
% S_obs = zeros(2000,9);
% parfor i = 1:2000 for parallel computing
% for i = 1:2000
%     [Pv, t] = sim_turin_matrix(N, Bw, Ns, theta_true);
%      For GPU acceleration
%     [Pv, t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_true); 
%      S_obs(i,:) = create_statistics(Pv, t);
% end

% Calculate mean and covariance of summaries on observed data
mu_S_obs = mean(S_obs);     % Mean of the summary statistics
Sigma_S_obs = cov(S_obs);     % Covariance of summary statistics

%% --- Initial max/min conditions for parameters (prior distribution) -----------------------------
load('Prior_min_max_values.mat')

%% --- ABC PMC algorithm ---------------------------------------------------------------------
% Number of summary statistics sets to generate
sumstat_iter = 2000;
% Extract this amount of parameter sets closest to observed data
nbr_extract = 100;
% Acceptance rate
Eps_percent = (nbr_extract/sumstat_iter)*100;
% Probability factor for generating population pool
prob_factor = 1e4;
% Number of ABC num_iter
num_iter = 4;

params_T = zeros(num_iter,nbr_extract);
params_G0 = zeros(num_iter,nbr_extract);
params_lambda = zeros(num_iter,nbr_extract);
params_sigma_N = zeros(num_iter,nbr_extract);

disp('BSL algorithm computing, please wait... ')

%% Iteration 1 - rejection step on uniform prior
out = zeros(5,sumstat_iter);
d = zeros(sumstat_iter,1);

% STEP 1: Sample parameter from predefined prior distribution (uniform):
% T (Reverberation time):
param_T = prior(1,1) + (prior(1,2) - prior(1,1)).*rand(sumstat_iter,1); % draw random numbers from prior
% G0 (Reverberation gain)
param_G0 = prior(2,1) + (prior(2,2) - prior(2,1)).*rand(sumstat_iter,1); % draw random numbers from prior
% lambda ()
param_lambda = prior(3,1) + (prior(3,2) - prior(3,1)).*rand(sumstat_iter,1); % draw random numbers from prior
% sigma_N (Variance noise floor)
param_sigma_N = prior(4,1) + (prior(4,2) - prior(4,1)).*rand(sumstat_iter,1); % draw random numbers from prior

% For parallel computing: parfor i = 1:sumstat_iter
for i = 1:sumstat_iter
    theta_curr = [param_T(i) param_G0(i) param_lambda(i) param_sigma_N(i)];
    
    % STEP 2: Simulate data using Turing model, based on parameters from STEP 1 and create statistics
    [Pv, t] = sim_turin_matrix(N, Bw, Ns, theta_curr);
    % For GPU acceleration
    % [Pv, t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_curr);
    
    % STEP 3: calculate the difference between observed and simulated summary statistics
    % Mahalanobis distance see formular in document.
    d(i) = (S_simulated - mu_S_obs)/Sigma_S_obs * (S_simulated - mu_S_obs)';
    
    % Row 1 of the out vector contains the distance
    % the rest of the rows contains the corresponding parameters
    % used for generating that specific distance.
    out(:,i) =  [d(i);...
        theta_curr'];
    disp(i);
end
% Sort the "out" matrix so that the lowest euclidean distance is at the
% (1,1) matrix position and highest distance is at (max,1)
out = sortrows(out',1)';
%
% Following vectors holds ALL extracted parameter values
% that was within euclidean distance
params_T(1,:)       = out(2,1:nbr_extract);
params_G0(1,:)      = out(3,1:nbr_extract);
params_lambda(1,:)  = out(4,1:nbr_extract);
params_sigma_N(1,:) = out(5,1:nbr_extract);
accepted_params = [params_T; params_G0; params_lambda; params_sigma_N];

%% sequential ABC (PMC-part)

% Calculate first weights and covariance
weights_T = ones(1,nbr_extract)./nbr_extract;
weights_G0 = ones(1,nbr_extract)./nbr_extract;
weights_lambda = ones(1,nbr_extract)./nbr_extract;
weights_sigma_N = ones(1,nbr_extract)./nbr_extract;
weights = [weights_T; weights_G0; weights_lambda; weights_sigma_N];
var_T = var(params_T(1,:));
var_G0 = var(params_G0(1,:));
var_lambda = var(params_lambda(1,:));
var_sigma_N = var(params_sigma_N(1,:));
covariance = diag([var_T var_G0 var_lambda var_sigma_N]);

% Choose theta from accepted parameters of last iteration with propbability based on wieghts
index_T = randsample((1:nbr_extract),sumstat_iter,true,weights(1,:));
index_G0 = randsample((1:nbr_extract),sumstat_iter,true,weights(2,:));
index_lambda = randsample((1:nbr_extract),sumstat_iter,true,weights(3,:));
index_sigma_N = randsample((1:nbr_extract),sumstat_iter,true,weights(4,:));

theta_prop = [params_T(1,index_T); params_G0(1,index_G0); params_lambda(1,index_lambda); params_sigma_N(1,index_sigma_N)];

for a = 2:num_iter%num_iter
    out = zeros(5,sumstat_iter);
    d = zeros(sumstat_iter,1);
    % For parallel computing: parfor i = 1:sumstat_iter
    for i = 1:sumstat_iter
        % Perturb theta
        theta_curr = mvnrnd(theta_prop(:,i),covariance);
        while(check_params(theta_curr,prior)==2)
            theta_curr = mvnrnd(theta_prop(:,i),covariance);
        end
        theta_curr(3) = round(theta_curr(3));
        %% STEP 2: Simulate data using Turing model, based on parameters from STEP 1 and create statistics
        [Pv, t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_curr);
        S_simulated = create_statistics(Pv, t);
        %% STEP 3: calculate the difference between observed and simulated summary statistics
        % Mahalanobis distance
        d(i) = (S_simulated - mu_S_obs)/Sigma_S_obs * (S_simulated - mu_S_obs)';
        
        % Row 1 of the out vector contains the distance
        % the rest of the rows contains the corresponding parameters
        % used for generating that specific distance.
        out(:,i) =  [d(i);...
            theta_curr'];
        disp(i);
    end
    % Sort the "out" matrix so that the lowest distance is at the
    % (1,1) matrix position and highest distance is at (1,end)
    out = sortrows(out',1)';
    
    % Following vectors holds ALL extracted parameter values
    % that was within euclidean distance
    params_T(a,:)       = out(2,1:nbr_extract);
    params_G0(a,:)      = out(3,1:nbr_extract);
    params_lambda(a,:)  = out(4,1:nbr_extract);
    params_sigma_N(a,:) = out(5,1:nbr_extract);
    
    accepted_params = [params_T(a,:); params_G0(a,:); params_lambda(a,:); params_sigma_N(a,:)];
    
    % Compute weights for next iteration
    old_weights = weights; % Save current weights
    for l = 1:size(weights,1)
        for k = 1:size(weights,2)
            weights(l,k) = 1 /(sum(pdf('normal',accepted_params(l,:), accepted_params(l,k), sqrt(covariance(l,l)).* old_weights(l,:))));
        end
    end
    weights = weights./sum(weights,2);
    % Find oovariance for next iteration
    covariance = 2*diag(diag(cov(out(2:5,1:nbr_extract)')));
    % Generate vector with each entry proportional weight for sampling
    probs = round(weights*prob_factor);
    % Generate population pool for sampling of new parameters
    big_T = [];
    big_G0 = [];
    big_lambda = [];
    big_sigma_N = [];
    for j = 1:nbr_extract
        big_T = [big_T repelem(accepted_params(1,j), probs(1,j))];
        big_G0 = [big_G0 repelem(accepted_params(2,j), probs(2,j))];
        big_lambda = [big_lambda repelem(accepted_params(3,j), probs(3,j))];
        big_sigma_N = [big_sigma_N repelem(accepted_params(4,j), probs(4,j))];
    end
    
    % Sample new para,eter values from population pool
    theta_prop = [datasample(big_T, sumstat_iter);...
        datasample(big_G0, sumstat_iter);...
        datasample(big_lambda, sumstat_iter);...
        datasample(big_sigma_N, sumstat_iter)];
    
    disp(a); % display iteration number
end
disp('PMC-ABC algorithm finished.')


%% Plots
thetas = [params_T; params_G0; params_lambda; params_sigma_N];

tt = tiledlayout(num_iter,4,'TileSpacing','Compact','Padding','Compact');
linesize = 2;
truecolor = '#32CD32';

for i = 1:num_iter
    [f_T, xi_T] = ksdensity(params_T(i,:));
    [f_G0, xi_G0] = ksdensity(pow2db(params_G0(i,:)));
    [f_lambda, xi_lambda] = ksdensity(params_lambda(i,:));
    [f_sigma_N, xi_sigma_N] = ksdensity(params_sigma_N(i,:));
    
    nexttile
    area(xi_T*1e9,f_T,'FaceColor','#bbbbbb','LineStyle','None');
    MMSE_T = mean(params_T(i,:));
    xline(MMSE_T*1e9,'r','LineWidth',linesize)
    xline(theta_true(1)*1e9,'--','Color',truecolor,'LineWidth',linesize)
    subtitle('T \times 10^{-9}')
    xlim([prior(1,1)*1e9 prior(1,2)*1e9])
    %set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    if i == 1
        legend( "Approx. posterior",'MMSEE',"True value")
    end
    ylabel(['Iteration: ' num2str(i)])
    nexttile
    area(xi_G0,f_G0,'FaceColor','#bbbbbb','LineStyle','None')
    xline(pow2db(theta_true(2)),'--','Color',truecolor,'LineWidth',linesize)
    MMSE_G0 = mean(params_G0(i,:));
    xline(pow2db(MMSE_G0),'r','LineWidth',linesize)
    subtitle('G_0 [dB]')
    xlim([pow2db(prior(2,1)) pow2db(prior(2,2))])
    %set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    nexttile
    area(xi_lambda*1e-9,f_lambda,'FaceColor','#bbbbbb','LineStyle','None')
    xline(theta_true(3)*1e-9,'--','Color',truecolor,'LineWidth',linesize)
    MMSE_lambda = mean(params_lambda(i,:));
    xline(MMSE_lambda*1e-9,'r','LineWidth',linesize)
    subtitle('\lambda \times 10^{9}')
    xlim([prior(3,1)*1e-9 prior(3,2)*1e-9])
    %set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    nexttile
    area(xi_sigma_N.^2*1e9,f_sigma_N,'FaceColor','#bbbbbb','LineStyle','None')
    xline(theta_true(4)^2*1e9,'--','Color',truecolor,'LineWidth',linesize)
    MMSE_sigma_N = mean(params_sigma_N(i,:));
    xline(MMSE_sigma_N^2*1e9,'r','LineWidth',linesize)
    subtitle('\sigma_N^2 \times 10^{-9}')
    xlim([prior(4,1).^2*1e9 prior(4,2).^2*1e9])
    %set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
end