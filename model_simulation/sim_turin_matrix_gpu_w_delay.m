function [P_y, t] = sim_turin_matrix_gpu_w_delay(N, B, Ns, Theta, delay)
    %  N       = Number of simulation iterations (The number of complete turin realizations)
    %  Bw      = Bandwidth 
    %  Ns      = Number of simulation iterations (The number of Hk realizations used) 
    %  T       = Reverberation time 
    %  G0      = Reverberation gain (power at delay zero)
    %  lambda  = Arrival rate of multipath components
    %  sigma_N = Noise variance
    
    T = Theta(1);
    G0 = Theta(2);
    lambda = Theta(3);
    sigma_N = Theta(4); 
    
    deltaf = B/(Ns-1); % Frequency seperation
    tmax = 1/deltaf; % Maximum delay, found from the bandwidth via frequency seperation
    t = (0:Ns-1)'./(deltaf*Ns); % Generate timestamps, in seconds
    %% Simulate model
    % We then simulate the transfer function, H_k, of the channel, using the
    % Turin model. 
    H_k = gpuArray(zeros(Ns,N)); % buffer for generated channel response data                                             
    % We run the simulation N times, creating new data sets for each
    % realization. 
    
    delay = tmax-delay; 
    f = B*(-Ns/2:(Ns/2-1))/Ns;
    shift = exp(-(1j*2*pi*f*delay))';
    
    parfor n = 1:N
        lmax = poissrnd(tmax*lambda);   % Number of multipath components, created from the Poisson distribution.
        tau = gpuArray(rand(lmax,1)*tmax);    % time-delays, drawn uniformly between 0 and the maximum delay.  
        % For every multipath component a complex gain is generated, based on a
        % sigma generated from a corresponding delay time value. 
        % The complex number is generated in cartesian form by drawing the real
        % and the imaginary part seperately from a normal distribution. 
        sigma_alpha = sqrt(G0*exp(-(tau*(1/T)) ) / lambda);% Calculate variance using eq 13 and Ph = Lambda*sigma^2 from Ayush paper.
        % The complex valued alpha is created by combining the real and
        % imaginary parts. 
        alpha = sigma_alpha * 1/sqrt(2) .* (randn(lmax,1) +  1j*randn(lmax,1));  
        % For every frequency index, k, the contribution from every multipath
        % component is added to the transfer function. 
        k = gpuArray((1:Ns));
        H_k(:,n) = (exp(-1j*2*pi*deltaf*k.*tau).' * alpha) .* shift;
    end
    % Generate noise vector
    noise = gpuArray(sigma_N^2 * (1/sqrt(2))* (randn(Ns,N) + 1j*randn(Ns,N)));
    % Add noise to transfer function with complex normal dist
    Y_k = (H_k + noise);
    % Power delay profile
    P_h = abs(ifft(Y_k,[],1)).^2;
    % We use P_Y = E_s * P_h + noise (Noise is already included in simulation)    
    P_y = P_h*B;
end