function val = check_params(params, prior)
% Function for checking validity of proposed parameter
    if        (params(1) <= prior(1,2) && params(1) >= prior(1,1) ...
            && params(2) <= prior(2,2) && params(2) >= prior(2,1) ...
            && params(3) <= prior(3,2) && params(3) >= prior(3,1) ...
            && params(4) <= prior(4,2) && params(4) >= prior(4,1))
        val = 1; 
    else
        val = 2;
    end
end