%Written by Brady Sivey

%The purpose of this function is to compute the necessary Jacobians
%to be used in the EKF.
function [Hk, Ac, Gamc] = EKFJacobians(xhatm, t, rE, wE, mu, stationvisindex)
%inputs
%1. xhatm - predicted state from dynamics
%2. t - time
%3. rE - Radius of Earth
%4. wE - Earth rotation rate
%5. mu - Gravitational Parameter
%6. stationvisindex - Vector of visible station numbers

%outputs
%1. Hk - Linearized measurement matrix for each station
%2. Ac - df/dx Jacobian
%3. Gamc - df/dw Jacobian

X = xhatm(1);
Xdot = xhatm(2);
Y = xhatm(3);
Ydot = xhatm(4);

N = length(stationvisindex);

Hk = zeros(3*N, 4);

%for loop to construct Hk
for i = 1:N
s = stationvisindex(i);

%givens
thetai0 = (s - 1) * pi/6;
Xs = rE*cos(wE*t + thetai0);
Ys = rE*sin(wE*t + thetai0);
Xsdot = -rE*wE*sin(wE*t + thetai0);
Ysdot = rE*wE*cos(wE*t + thetai0);

%These equations come from my derivations
dX = X - Xs;
dXdot = Xdot - Xsdot;
dY = Y - Ys;
dYdot = Ydot - Ysdot;

rho = sqrt(dX^2 + dY^2);
A = dX*dXdot + dY*dYdot;

% 3x4 H matrix
H21 = (rho*dXdot - A*(dX/rho))/(rho^2);
H23 = (rho*dYdot - A*(dY/rho))/(rho^2);
Hi = [dX/rho 0 dY/rho 0;
    H21 dX/rho H23 dY/rho;
    -dY/(rho^2) 0 dX/(rho^2) 0];

idx = 3*(i-1) + (1:3);
Hk(idx,:) = Hi;

end

%constructing Ac and Gamc matrices

%These come from my derivations
Ac21 = -mu * (Y^2 - 2*X^2)/(X^2 + Y^2)^(5/2);
Ac23 = mu * (3*X*Y)/(X^2 + Y^2)^(5/2);
Ac41 = mu * (3*X*Y)/(X^2 + Y^2)^(5/2);
Ac43 = -mu * (X^2 - 2*Y^2)/(X^2 + Y^2)^(5/2);

Ac = [0 1 0 0;
    Ac21 0 Ac23 0
    0 0 0 1;
    Ac41 0 Ac43 0];

Gamc = [0 0;
    1 0;
    0 0;
    0 1];

end