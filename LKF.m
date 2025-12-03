function [x_full,y_full,Ppkp1] = LKF(F,G,H,Q,R,Omega,delta_x0,Pp0,x_nom,u_nom,u,y_nom,y)
    % Initialize
    N = size(F,3);
    delta_xpkp1 = zeros(4,N);
    delta_xpkp1(:,1) = delta_x0;
    delta_ykp1 = cell(N,1);
    y_full = cell(N,1);
    Ppkp1(:,:,1) = Pp0;
    delta_uk = u - u_nom;
   
    % LKF Algorithm
    for k = 1:N-1
        % Time update/prediction step
        delta_xmkp1 = F(:,:,k)*delta_xpkp1(:,k) + G*delta_uk(:,k);
        Pmkp1 = F(:,:,k)*Ppkp1(:,:,k)*F(:,:,k)' + Omega*Q*Omega';

        % Measurement update/correction step
        K = height(H{k+1})/3;
        Rk = kron(eye(K),R);
        Kkp1 = Pmkp1*H{k+1}'*inv(H{k+1}*Pmkp1*H{k+1}'+Rk);
        delta_ykp1{k+1} = y{k+1} - y_nom{k+1};
        Ppkp1(:,:,k+1) = (eye(4)-Kkp1*H{k+1})*Pmkp1;
        delta_xpkp1(:,k+1) = delta_xmkp1+Kkp1*(delta_ykp1{k+1}-H{k+1}*delta_xmkp1);


        % Add Outputs
        y_full{k+1} = y_nom{k+1} + delta_ykp1{k+1};
    end

    % Add to nominal state
    x_full = x_nom' + delta_xpkp1;
end

