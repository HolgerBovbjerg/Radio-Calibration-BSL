%% clear
clear 

%% Initial choices made on model. 
% We have chosen some values based on "Estimator for Stochastic Channel Model without
% Multipath Extraction using Temporal Moments" by Ayush Bharti et al. 

N = 200;% Number of realizations to generate and average over.
        % This control the variance of the estimate, i.e. more realizations
        % less variance

B = 4e9; % Bandwidth of signal: 4 [GHz]
Ns = 801; % Number of sample points in each data set
T = 7.8e-9; % Reverberation time: 7.8 [ns]
G0 = db2pow(-83.9); % Reverberation gain -83.9 [dB] converted from dB to linear scale
lambda = 1e9; % arrival rate lambda 1 [GHz]
sigma_N = sqrt(0.28e-9); % Noise standard deviation
% Collect parameters in theta
theta = [T G0 lambda sigma_N];

% Generate realizations of Turin model
[P_y, t] = sim_turin_matrix(N,B,Ns,theta); % N realisations of Turin model power delay profile
% For gpu acceleration
% [P_y, t] = sim_turin_matrix_gpu(N,B,Ns,theta); % N realisations of Turin model power delay profile

% We use the the reverberation model as theoretical comparison 
P_h_theoretical = G0*exp(-(t/T));

% We use P_Y = E_s * P_h + noise (Noise is already included in simulation)
P_y_simulated = mean(P_y,2); 
P_y_theoretical = P_h_theoretical + sigma_N^2/Ns; % theoretical does not need bandwidth scaling


%% Generation of plots showing the power spectrum. 
figure
plot(t*1e9,pow2db(P_y_theoretical), 'DisplayName', "P_y theoretical")
hold on
plot(t*1e9,pow2db(P_y_simulated), 'DisplayName', "P_y simulated")
xlim([0 200])
ylim([-130 -80])
xlabel("Time [ns]")
ylabel("Power [dB]")
lgd = legend;