clc; clear all;

%% Simulate Non-Linear ODEs and Measured Outputs

% Constants
mu = 398600; % km^3/s^2
r0 = 6678; % km
Ydot0 = r0 * sqrt(mu/r0^3); % km/s
RE = 6378; % km
wE = (2*pi)/86400; % rad/s
w = sqrt(mu/r0^3);

% ICs
nom_var0 = [r0 0 0 Ydot0]';
delta_x0 = [0 0.075 0 -0.021]';
pert_var0 = nom_var0 + delta_x0;

% Simulation time
delta_t = 10; %s
tspan = 0:delta_t:14000; % s

%% CT Dynamics
Abar = @(t) [0 1 0 0
        -w^2*(1-3*cos(w*t)^2) 0 (3/2)*w^2*sin(2*w*t) 0
        0 0 0 1
        (3/2)*w^2*sin(2*w*t) 0 -w^2*(1-3*sin(w*t)^2) 0];

Bbar = [0 0
        1 0
        0 0
        0 1];

thetai0 = @(i) (i - 1) * pi/6;
Xi = @(t,i) RE * cos(wE*t + thetai0(i));
Yi = @(t,i) RE * sin(wE*t + thetai0(i));
thetai = @(t,i) atan2(Yi(t,i),Xi(t,i));
Xidot = @(t,i) -wE * RE * sin(wE*t + thetai0(i));
Yidot = @(t,i) wE * RE * cos(wE*t + thetai0(i));

delta_Xnom = @(t,i) r0*cos(w*t) - Xi(t,i);
delta_Ynom = @(t,i) r0*sin(w*t) - Yi(t,i);
delta_Xdotnom = @(t,i) -w*r0*sin(w*t) - Xidot(t,i);
delta_Ydotnom = @(t,i) w*r0*cos(w*t) - Yidot(t,i);
Anom = @(t,i) delta_Xnom(t,i)*delta_Xdotnom(t,i)...
    + delta_Ynom(t,i)*delta_Ydotnom(t,i);
rhonom = @(t,i) sqrt(delta_Xnom(t,i)^2 + delta_Ynom(t,i)^2);

Cbar = @(t,i) [delta_Xnom(t,i)./rhonom(t,i), 0, delta_Ynom(t,i)./rhonom(t,i)...
    , 0; % row 1
(rhonom(t,i)*delta_Xdotnom(t,i) - Anom(t,i)*(delta_Xnom(t,i)/rhonom(t,i)))...
/ rhonom(t,i)^2,... row 2 col 1
delta_Xnom(t,i)/rhonom(t,i),... row 2 col 2
(rhonom(t,i)*delta_Ydotnom(t,i) - Anom(t,i)*(delta_Ynom(t,i)/rhonom(t,i)))...
/ rhonom(t,i)^2,... row 2 col 3
delta_Ynom(t,i)/rhonom(t,i); % row 2 col 4
-delta_Ynom(t,i)/rhonom(t,i)^2, 0, delta_Xnom(t,i)/rhonom(t,i)^2, 0]; % row 3

%% Load in Data Logs and Define Inputs
% Load in Data Logs and Q/R
data = load('orbitdeterm_finalproj_KFdata.mat');
Q = data.Qtrue;
R = data.Rtrue;
tvec_datalog = data.tvec;
ydatalog = data.ydata;

% ode45 call
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode45(@(tspan,var) OrbitEOM(tspan,var,mu),tspan,pert_var0,options);
[~,state_nom] = ode45(@(tspan,var0) OrbitEOM(tspan,var0,mu),tspan,nom_var0,options);

% Get outputs
[output_var,station_vis] = MeasuredOutput(t,state,RE,wE,true);
[output_var_nom,~] = MeasuredOutput(t,state_nom,RE,wE,false);

% Convert into cell array - filter nominal measurements to match NL idxs
N = length(t);
ycell_nom = cell(N,1);
ycell_pert = cell(N,1);
station_vis_cell = cell(N,1);
for k = 1:N-1
    yrowsKeep = ~any(isnan(output_var(:,k+1)),2); % rows that contain NO NaNs
    stationrowsKeep = ~any(isnan(station_vis(:,k+1)),2); % rows that contain NO NaNs
    y1 = output_var_nom(yrowsKeep,k+1);
    y2 = output_var(yrowsKeep,k+1);
    station = station_vis(stationrowsKeep,k+1);
    ycell_nom{k+1} = y1;
    ycell_pert{k+1} = y2;
    station_vis_cell{k+1} = station;
end

y_pert = ycell_pert;
y_nom = ycell_nom;
station_vis = station_vis_cell;

% Concatenate Outputs/Stations
NL_outputs = [y_nom station_vis];

Dbar = zeros(3,2);

n_sims = 10;
% epsilon_x = zeros(size(x,2),n_sims);
for sim_num = 1:n_sims
    %% CT -> DT and Simulate Perturbation Dynamics
    % Initialize
    delta_xk = zeros(4, N);
    Pp0 = diag([10,0.1,10,0.1]);
    delta_xk(:,1) = chol(Pp0,'lower') * randn(4,1);
    delta_yk = cell(N,1);
    Fk = zeros(4,4,N);
    Hk = cell(N,1); 
    yk = cell(N,1);
    
    % Process Noise Matrix
    Gamma = Bbar; % w1,2 has same mapping as u1,2
    Omegabar = delta_t.*Gamma;
    Sw = chol(Q,'lower');
    for k = 1:N-1
        % Build H at time t_k
        num_vis = station_vis{k+1};
        H = zeros(3*length(num_vis),4);
        for i = 1:length(num_vis)
            H(3*i-3+(1:3),:) = Cbar(tspan(k+1),num_vis(i));
        end
        Hk{k+1} = H;
    
        % Propagate perturbation to next step (except at the last time)
        if k < N
            Fk(:,:,k) = eye(4) + delta_t.*Abar(tspan(k));
            G = delta_t.*Bbar;                  
            delta_xk(:,k+1) = Fk(:,:,k)*delta_xk(:,k) + Omegabar*Sw*randn(2,1);
        end
    
        % Linearized measurement at time t_k
        delta_yk{k+1} = Hk{k+1}*delta_xk(:,k+1);
        yk{k+1} = y_nom{k+1} + delta_yk{k+1};
    end
    
    % Reconstruct the linearized full state
    state_lin = state_nom + delta_xk.';
    
    % Concatanate Outputs/Stations
    L_outputs = [yk,station_vis];
    
    % Inputs
    u_nom = zeros(2,length(tspan));
    u = zeros(2,length(tspan));
    
    %% Noise and Covariance
    % Generate Noisy Measurements;
    y_pert_noise = cell(N,1);
    for k = 1:N-1
        K = length(y_pert{k+1});
        Sv = chol(kron(eye(K/3),R),'lower');
        qk = randn(K,1);
        y_pert_noise{k+1} = y_pert{k+1} + Sv*qk;
    end
    
    %% LKF
    lkf_delta_x0 = zeros(4,1);
    [x_LKF,y_LKF,Ppkp1_LKF] = LKF...
        (Fk,G,Hk,Q,R,Omegabar,lkf_delta_x0,Pp0,state,u_nom,u,y_nom,y_pert_noise);
    
    LKF_outputs = [y_LKF station_vis];

    epsilon_i_x = zeros(size(state_lin,2),1); % kx1
    for k = 1:N 
        x_k = state_lin(k,:)';
        xhatp_k = x_LKF(:,k);
        Pp_k = Ppkp1_LKF(:,:,k);
        e_xk = x_k - xhatp_k;
        epsilon_i_xk = e_xk'*Pp_k^-1*e_xk; % 1x1
        epsilon_i_x(k,1) = epsilon_i_xk;
    end
    epsilon_x(:,sim_num) = epsilon_i_x;

end

epsilonbar_x = sum(epsilon_x,2)/n_sims;
figure;
scatter(1:1401, epsilonbar_x, 'filled'); grid on;
xlabel('Time Step k');
ylabel('$\bar{\epsilon}_x$', 'Interpreter', 'latex');
title('LKF NEES Test Results');