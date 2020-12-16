function [P_y, t] = sim_turin_matrix(N, Bw, Ns, Theta)
    %% INPUTS:
    %  N       = Number of simulation iterations (The number of complete turin realizations)
    %  Bw      = Bandwidth [Hz]
    %  Ns      = Number of samples for each realization
    %  T       = Reverberation time [s]
    %  G0      = Reverberation gain (power at delay zero)
    %  lambda  = Mean arrival rate of multipath components per second [Hz]
    %  sigma_N = Noise standard deviation
   
    %% OUTPUTS:
    % P_Y       = Power delay spectrum of model realizations
    % t         = Time vector
 
    % Model parameters input
    T = Theta(1);
    G0 = Theta(2);
    lambda = Theta(3);
    sigma_N = Theta(4);

    % Calculate signal values 
    deltaf = Bw/(Ns-1); % Frequency seperation
    tmax = 1/deltaf; % Maximum delay, found from the bandwidth via frequency seperation
    t = (0:Ns-1)'./(deltaf*Ns); % Generate timestamps, in seconds
    
    %% Simulate model
    % Simulate the transfer function, H_k, of the channel, using the Turin model. 
    H_k = zeros(Ns,N); % buffer for generated channel response data                                             
    % Run the simulation N times, creating N realizations of the model. 
    % For parallel computing: parfor n = 1:N
    for n = 1:N    
        lmax = poissrnd(tmax*lambda); % Number of multipath components, created from the Poisson distribution.
        tau = rand(lmax,1)*tmax; % time-delays, drawn uniformly between 0 and the maximum delay.  
        % For every multipath component a complex gain is generated, based on a sigma generated from a corresponding delay time value. 
        % The complex number is generated in cartesian form by drawing the real and the imaginary part seperately from a normal distribution. 
        sigma_alpha = sqrt(G0*exp(-(tau*(1/T)) ) / lambda);% Calculate variance using eq 13 and Ph = Lambda*sigma^2 from Ayush paper.
        % The complex valued alpha is created by combining the real and imaginary parts.
        alpha = sigma_alpha * 1/sqrt(2) .* (randn(lmax,1) +  1j*randn(lmax,1));
        % For every frequency index, k, the contribution from every multipath component is added to the transfer function. 
        k = (1:Ns);
        H_k(:,n) = (exp(-1j*2*pi*deltaf*k.*tau).' * alpha); % (Hk(:,n) - the nth column of matrix Hk)
    end
    
    % Generate noise vector from complex normal distribution with mean 0 and variance sigma_N^2
    noise = sigma_N^2 * (1/sqrt(2))* (randn(Ns,N) + 1j*randn(Ns,N));

    % Add noise to transfer function
    Y_k = H_k + noise;


    % Power delay profile:
    % Inverse furier transform of each column of Hk and then the absolute value
    % of this, elementwise squared. P_y -> A definition from Ayush paper.
    % P_y is now the Power delay profile in the time domain
    P_y = abs(ifft(Y_k,[],1)).^2;
    % We use P_Y = E_s * P_h + noise (Noise is already included in simulation)
    % E_s = energy of signal 
    P_y = Bw*P_y; % Bandwidth scaling
    
end