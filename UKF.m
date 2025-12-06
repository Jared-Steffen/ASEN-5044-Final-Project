function [outputArg1,outputArg2] = UKF(t,xp0,P0)
    % xp0: initial state condition
    % P0: inital covariance matrix

    % Initialize/Preallocate
    xpkp1 = zeros(length(t),length(P0));
    Ppkp1 = zeros(size(P0),length(t));

end

