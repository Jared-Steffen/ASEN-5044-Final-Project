function [outputArg1,outputArg2] = LKF(F,G,H,M,Q,R,delta_x0,delta_u0,P0,x_nom,u_nom)

    % Initialize
    delta_xk = zeros(4, N);
    delta_xk(:,1) = delta_x0;
    delta_uk = zeros(2, N);
    delta_uk(:,1) = delta_u0;

   

    for k = 1:size(F,3)

        % Time update/prediction step
        delta_xkp1m = F(:,:,k)*delta_xk(:,k) + G*delta_uk;

end

