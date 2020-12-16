function [covariance, thetacurr] = find_cov_prior(prior)

T_var = 1/12 * (prior(1,2) - prior(1,1))^2;
G0_var = 1/12 * (prior(2,2) - prior(2,1))^2;
lambda_var = 1/12 * (prior(3,2) - prior(3,1))^2;
sigmaN_var = 1/12 * (prior(4,2) - prior(4,1))^2;

T_prior = prior(1,1) + (prior(1,2) - prior(1,1)).*rand;
G0_prior = prior(2,1) + (prior(2,2) - prior(2,1)).*rand;
lambda_prior = prior(3,1) + (prior(3,2) - prior(3,1)).*rand;
sigmaN_prior = prior(4,1) + (prior(4,2) - prior(4,1)).*rand;

thetacurr = [T_prior G0_prior lambda_prior sigmaN_prior]';
covariance = diag([T_var G0_var lambda_var sigmaN_var]);

end
