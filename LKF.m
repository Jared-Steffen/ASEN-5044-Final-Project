function [x_full,y_full,Ppkp1,innov,delta_xpkp1,Skp1] = LKF...
    (F,G,H,Q,R,Omega,delta_x0,Pp0,x_nom,u_nom,u,y_nom,y_noisy)

    % Preallocate/initialize
    N = size(F,3);
    delta_xpkp1 = zeros(4,N);
    delta_xpkp1(:,1) = delta_x0;
    innov = cell(N,1);
    y_full = cell(N,1);
    Ppkp1(:,:,1) = Pp0;
    delta_uk = u - u_nom;
    Skp1 = cell(N,1);
    Skp1{1} = zeros(length(y_noisy{1}), length(y_noisy{1}));
   
    % LKF Algorithm
    for k = 1:N-1
        % Time update/prediction step
        delta_xmkp1 = F(:,:,k)*delta_xpkp1(:,k) + G*delta_uk(:,k);
        Pmkp1(:,:,k) = F(:,:,k)*Ppkp1(:,:,k)*F(:,:,k)' + Omega*Q*Omega';

        % Measurement update/correction step
        K = height(H{k+1})/3;
        Rk = kron(eye(K),R);
        Skp1{k+1} = (H{k+1}*Pmkp1(:,:,k)*H{k+1}'+Rk);
        Kkp1 = Pmkp1(:,:,k)*H{k+1}'/Skp1{k+1};
        delta_ykp1 = y_noisy{k+1} - y_nom{k+1};
        Ppkp1(:,:,k+1) = (eye(4)-Kkp1*H{k+1})*Pmkp1(:,:,k);
        delta_xpkp1(:,k+1) = delta_xmkp1+Kkp1*(delta_ykp1-H{k+1}*delta_xmkp1);

        % Innovation
        innov{k+1} = delta_ykp1-H{k+1}*delta_xmkp1;

        % Add Outputs
        y_full{k+1} = y_nom{k+1} + H{k+1}*delta_xmkp1;
    end

    % Add to nominal state
    x_full = x_nom + delta_xpkp1';
end