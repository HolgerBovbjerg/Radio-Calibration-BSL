%% ABC implementation - ABC REJECTION ALGORITHM:
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
%     [Pv, t] = sim_turin_matrix_gpu_gpu(N, Bw, Ns, theta_true);
%      S_obs(i,:) = create_statistics(Pv, t);
% end

% Calculate mean and covariance of summaries on observed data
mu_S_obs = mean(S_obs);     % Mean of the summary statistics 
Sigma_S_obs = cov(S_obs);     % Covariance of summary statistics

%% --- Initial max/min conditions for parameters (prior distribution) -----------------------------
load('Prior_min_max_values.mat')

%% --- ABC rejection algorithm ---------------------------------------------------------------------
% Number of summary statistics sets to generate  
sumstat_iter = 20;
% Extract this amount of parameter sets closest to observed data
nbr_extract = 5;
% Acceptance rate
Eps_percent = (nbr_extract/sumstat_iter)*100;
 

disp('ABC algorithm computing, please wait... ')

out = zeros(5,sumstat_iter);
d = zeros(sumstat_iter,1);

%% STEP 1: Sample parameter from predefined prior distribution (uniform):
% T (Reverberation time):
param_T = prior(1,1) + (prior(1,2) - prior(1,1)).*rand(sumstat_iter,1); % draw random numbers from prior
% G0 (Reverberation gain)
param_G0 = prior(2,1) + (prior(2,2) - prior(2,1)).*rand(sumstat_iter,1); % draw random numbers from prior 
% lambda ()
param_lambda = prior(3,1) + (prior(3,2) - prior(3,1)).*rand(sumstat_iter,1); % draw random numbers from prior 
% sigma_N (Variance noise floor)
param_sigma_N = prior(4,1) + (prior(4,2) - prior(4,1)).*rand(sumstat_iter,1); % draw random numbers from prior 


% for parallel computing: parfor i = 1:sumstat_iter
for i = 1:sumstat_iter
    theta_curr = [param_T(i) param_G0(i) param_lambda(i) param_sigma_N(i)]; % Current theta proposal
    
    %% STEP 2: Simulate data using Turing model, based on parameters from STEP 1 and create statistics
    [Pv, t] = sim_turin_matrix(N, Bw, Ns, theta_curr);
    % For GPU acceleration
    % [Pv, t] = sim_turin_matrix_gpu(N, Bw, Ns, theta_curr);
    S_simulated = create_statistics(Pv, t);
    %% STEP 3: calculate the distance between observed and simulated summary statistics 
    % Mahalanobis distance
    d(i) = (S_simulated - mu_S_obs)/Sigma_S_obs * (S_simulated - mu_S_obs)';

    % Row 1 of the out vector contains the distance 
    % the rest of the rows contains the corresponding parameters 
    % used for generating that specific distance.    
    out(:,i) =  [d(i);...
                theta_curr'];
    disp(i);
end

params_T       = out(2,1:nbr_extract);
params_G0      = out(3,1:nbr_extract);
params_lambda  = out(4,1:nbr_extract);
params_sigma_N = out(5,1:nbr_extract);

accepted_params = [params_T; params_G0; params_lambda; params_sigma_N];
disp('ABC algorithm finished.')

%% Plots
[f_T, xi_T] = ksdensity(accepted_params(1,:));
[f_G0, xi_G0] = ksdensity(pow2db(accepted_params(2,:)));
[f_lambda, xi_lambda] = ksdensity(accepted_params(3,:));
[f_sigma_N, xi_sigma_N] = ksdensity(accepted_params(4,:));

tt = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
linesize = 2;
truecolor = '#32CD32';

nexttile
area(xi_T*1e9,f_T,'FaceColor','#bbbbbb','LineStyle','None');
MMSE_T = mean(accepted_params(1,:));
xline(MMSE_T*1e9,'r','LineWidth',linesize)
xline(theta_true(1)*1e9,'--','Color',truecolor,'LineWidth',linesize)
subtitle('T \times 10^{-9}')
xlim([prior(1,1)*1e9 prior(1,2)*1e9])
%set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
legend( "Approx. posterior",'MMSE',"True value")

nexttile
area(xi_G0,f_G0,'FaceColor','#bbbbbb','LineStyle','None')
xline(pow2db(theta_true(2)),'--','Color',truecolor,'LineWidth',linesize)
MMSE_G0 = mean(accepted_params(2,:));
xline(pow2db(MMSE_G0),'r','LineWidth',linesize)
subtitle('G_0 [dB]')
xlim([pow2db(prior(2,1)) pow2db(prior(2,2))])
%set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

nexttile
area(xi_lambda*1e-9,f_lambda,'FaceColor','#bbbbbb','LineStyle','None')
xline(theta_true(3)*1e-9,'--','Color',truecolor,'LineWidth',linesize)
MMSE_lambda = mean(accepted_params(3,:));
xline(MMSE_lambda*1e-9,'r','LineWidth',linesize)
subtitle('\lambda \times 10^{9}')
xlim([prior(3,1)*1e-9 prior(3,2)*1e-9])
%set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

nexttile
area(xi_sigma_N.^2*1e9,f_sigma_N,'FaceColor','#bbbbbb','LineStyle','None')
xline(theta_true(4)^2*1e9,'--','Color',truecolor,'LineWidth',linesize)
MMSE_sigma_N = mean(accepted_params(4,:));
xline(MMSE_sigma_N^2*1e9,'r','LineWidth',linesize)
subtitle('\sigma_N^2 \times 10^{-9}')
xlim([prior(4,1).^2*1e9 prior(4,2).^2*1e9])
%set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])