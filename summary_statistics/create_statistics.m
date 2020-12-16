function S = create_statistics(Pv, t)
    % 0th, 1st and 2nd moments of Pv = |y|^2
    m0 = trapz(t,Pv)';         % 0th second order moment (t^0*Pv = Pv)
    m1 = trapz(t,t.^1.*Pv)';   % 1st second order moment
    m2 = trapz(t,t.^2.*Pv)';   % 2nd second order moment 
    
    % Create matrix of second order moements
    D2 = [m0 m1 m2];
    
    % Calculate mean and covariance 
    meanD2 = mean(D2);
    CovD2 = cov(D2);
    
    % Create summary statistic vector
    S = [meanD2(1) meanD2(2) meanD2(3)...
        CovD2(1,1) CovD2(2,2) CovD2(3,3)...
        CovD2(1,2) CovD2(1,3) CovD2(2,3)];

    % Check for negative statistics
    for i = 1:length(S)
        if S(i) < 0
            S(i) = 0;
        end
    end
    % Take natural logarithm
    S = log(S);
end
