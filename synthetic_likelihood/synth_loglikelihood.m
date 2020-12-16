function L = synth_loglikelihood(s, s_star)
    mu = mean(s_star); % Calculate mean of summary statistics from simulated data 
    Sigma = cov(s_star); % Calculate covariance of summary statistics from simulated data 
    %     L = -1/2*( (s - mu)*inv(Sigma))*(s - mu)' - 1/2*log(det(Sigma)); % Synthetic likelihood L
    L = -1/2*( ( (s - mu)/(Sigma)) ) * (s - mu)' - 1/2*log(det(Sigma)); % Synthetic likelihood L
end

