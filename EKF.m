%Written by Brady Sivey

%This function implements the Extended Kalman Filter as learned in ASEN
%5044
function [xhatp, Pp, yhat, Sv] = EKF(Q, R, ydata, tvec, mu, rE, wE, xhat0, P0, station_vis)
%inputs
%1. Q - Process noise covariance
%2. R - Measurement Noise Covariance
%3. ydata - Measurements
%4. tvec - Time vector
%5. mu - Gravitational Parameter
%6. rE - Radius of Earth
%7. wE - Earth rotation rate
%8. xhat0 - Initial state estimate
%9. P0 - Initial Covariance
%10. station_vis - Vector of visible station numbers

%outputs
%1. xhatp - state estimates
%2. Pp - Covariance
%3. yhat - Predicted Measurements

N = length(tvec);

%preallocating
xhatp = zeros(4, N);
Pp = zeros(4, 4, N);
yhat = cell(N,1);
Sv = cell(N,1);

%step 1
%initializing
xhatp(:, 1) = xhat0;
Pp(:, :, 1) = P0;

for k = 1:N-1
%step 2
tk = tvec(k);
tkn = tvec(k+1);
dt = tkn - tk;

%step 3
x0 = xhatp(:, k);
%deterministic nonlinear DT function eval
[~, xtraj] = ode45(@(t,x) OrbitEOM(t,x,mu), [tk tkn], x0);
xhatm = xtraj(end,:).';

vis_idx = station_vis{k+1};
[Hk, Ac, Gamc] = EKFJacobians(xhatm, tkn, rE, wE, mu, vis_idx);

Fk = eye(4) + dt * Ac;
Omegak = dt * Gamc;

%approx predicted covariance
Pm = Fk*Pp(:,:,k)*Fk.' + Omegak*Q*Omegak.';

%step 4
NaN_bool = true;
%deterministic nonlinear function eval
ttmp = [0; tkn];
statetmp = [xhatm.'; xhatm.']; 
[yfa, ~] = MeasuredOutput(ttmp, statetmp, rE, wE, NaN_bool);
yf = yfa(:,2);

m = length(vis_idx);
yhatf = zeros(3*m,1);
for j=1:m
s = vis_idx(j);
idxA = 3*(s-1) + (1:3);
idxV = 3*(j-1) + (1:3);
yhatf(idxV) = yf(idxA);
end
yhat{k+1} = yhatf;

%innovation
innov = ydata{k+1} - yhatf;

Rk = kron(eye(m), R);

%kalman gain
Sv{k+1} = (Hk * Pm * Hk.' + Rk);
K = Pm * Hk.' / Sv{k+1};

%updated total state estimate
xhatp(:,k+1) = xhatm + K*innov;

%updated covariance approximation
Pp(:,:,k+1) = (eye(4) - K*Hk) * Pm;

end