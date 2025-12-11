clc; clear; close all

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
% T = (2*pi)/(sqrt(mu/r0^3));
% tspan = 0:delta_t:T; % s
tspan = 0:delta_t:14000; % s

% ode45 call 
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[~,x_nom] = ode45(@(tspan,var0) OrbitEOM(tspan,var0,mu),tspan,nom_var0,options);

% ode45 call - pert 
[t,x_pert] = ode45(@(tspan,var) OrbitEOM(tspan,var,mu),tspan,pert_var0,options);

% Get outputs
[output_var,station_vis] = MeasuredOutput(t,x_pert,RE,wE,true);
[output_var_nom,~] = MeasuredOutput(t,x_nom,RE,wE,false);

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

% Concatanate Outputs/Stations
NL_outputs = [y_nom station_vis];

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

Dbar = zeros(3,2);

%% CT -> DT and Simulate Perturbation Dynamics
% Initialize
delta_xk = zeros(4, N);
delta_xk(:,1) = delta_x0;
delta_yk = cell(N,1);
Fk = zeros(4,4,N);
Hk = cell(N,1); 
yk = cell(N,1);

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
        delta_xk(:,k+1) = Fk(:,:,k)*delta_xk(:,k);
    end

    % Linearized measurement at time t_k
    delta_yk{k+1} = Hk{k+1}*delta_xk(:,k+1);
    yk{k+1} = y_nom{k+1} + delta_yk{k+1};
end

% Reconstruct the linearized full state
state_lin = x_nom + delta_xk.';

% Concatanate Outputs/Stations
L_outputs = [yk,station_vis];

%% Load in Data Logs and Define Inputs
% Load in Data Logs and Q/R
data = load('orbitdeterm_finalproj_KFdata.mat');
Qtrue = data.Qtrue;
R = data.Rtrue;
tvec_datalog = data.tvec;
ydatalog = data.ydata';
ydatalog(1) = cell(1,1);


N_log = length(ydatalog);
t_log = tvec_datalog(:);

ydatalog_mod        = cell(N_log,1);
station_vis_datalog = cell(N_log,1);

for k = 1:N_log
    current_y = cell2mat(ydatalog(k));

    if isempty(current_y)
        % No measurements at this time step
        ydatalog_mod{k} = zeros(0,1);
        station_vis_datalog{k} = [];
    else
        station_vis_datalog{k} = current_y(end,:)';
        ytmp = current_y(1:end-1,:);
        ydatalog_mod{k} = ytmp(:);
    end
end

y_nom_log = cell(N_log,1);
for k = 1:N_log
    stations = station_vis_datalog{k};

    if isempty(stations)
        y_nom_log{k} = zeros(0,1);
        continue;
    end

    yk_nom = zeros(3*length(stations),1);
    for i = 1:length(stations)
        s = stations(i); % station index
        rows = 3*s-2 : 3*s; % rows in output_var_nom for station s
        yk_nom(3*i-3+(1:3)) = output_var_nom(rows, k);
    end
    y_nom_log{k} = yk_nom;
end

% Inputs
u_nom = zeros(2,length(tspan));
u = zeros(2,length(tspan));

%% Tuning Knobs
Q_LKF = 10000*Qtrue;
Q_EKF = 1*Qtrue;
Q_UKF = 1*Qtrue;
alpha = 1e-4;
beta = 2;
kappa = 0;

%% Noise and Covariance
% Process Noise Matrix
Gamma = Bbar; % w1,2 has same mapping as u1,2
Omegabar = delta_t.*Gamma;

% Generate Noisy Measurements;
y_pert_noise = cell(N,1);
for k = 1:N-1
    K = length(y_pert{k+1});
    Sv = chol(kron(eye(K/3),R),'lower');
    qk = randn(K,1);
    y_pert_noise{k+1} = y_pert{k+1} + Sv*qk;
end

noisy_ouputs = [y_pert_noise station_vis];

% ode45 call - pert w/ noise
x_true_pert = pert_var0;
x_pert_noisy = zeros(N,4);
x_pert_noisy(1,:) = pert_var0';
for k = 1:N-1
    current_tspan = [tspan(k) tspan(k+1)];
    [~,x_pert_k] = ode45(@(current_tspan,x_true_pert)...
        OrbitEOM(current_tspan,x_true_pert,mu),current_tspan,x_true_pert,options);
    xtrue_kp1 = x_pert_k(end,:); 
    w_k = chol(Qtrue,"lower")*randn(2,1);
    x_true_pert = xtrue_kp1' + Omegabar*w_k;
    x_pert_noisy(k+1,:) = x_true_pert;
end

% Initialize Covariance
Pp0 = diag([100,1,100,1]);

%% LKF
[x_LKF,y_LKF,P_LKF,innov_LKF,delta_x_LKF,Sv_LKF] = LKF...
    (Fk,G,Hk,Q_LKF,R,Omegabar,delta_x0,Pp0,x_nom,u_nom,u,y_nom,y_pert_noise);

LKF_outputs = [y_LKF station_vis];
LKF_state_err = x_LKF-x_pert;

%% EKF
[x_EKF,P_EKF,y_EKF,innov_EKF,Sv_EKF] = EKF(Q_EKF,R,y_pert_noise,t,mu,RE,wE,nom_var0,Pp0,station_vis);

EKF_outputs = [y_EKF station_vis];
EKF_state_err = x_EKF'-x_pert;

%% UKF
[x_UKF,P_UKF,y_UKF,innov_UKF,Sv_UKF] = UKF(t,mu,RE,wE,nom_var0,Pp0,Q_UKF,...
    R,Omegabar,alpha,beta,kappa,y_pert_noise,station_vis);

UKF_outputs = [y_UKF station_vis];
UKF_state_err = x_UKF'-x_pert;

%% Plots
% Dynamics Labels
Full_Dynamics_Labels = {'$X$ [km]','$\dot{X}$ [km/s]','$Y$ [km]',...
    '$\dot{Y}$ [km/s]'};
Perturbation_Dynamics_Labels = {'$\delta X$ [km]',...
    '$\delta\dot{X}$ [km/s]','$\delta Y$ [km]','$\delta\dot{Y}$ [km/s]'};
Error_Dynamics_Labels = {'$X$ Error [km]','$\dot{X}$ Error [km/s]',...
    '$Y$ Error [km]','$\dot{Y}$ Error[km/s]'};

% NL System
Plot_Dynamics(t,x_pert,Full_Dynamics_Labels,'Nonlinear Dynamics')
Plot_Outputs(t,NL_outputs,'Nonlinear Model Outputs')

% Linearized System
Plot_Dynamics(t,delta_xk',Perturbation_Dynamics_Labels,...
    'Linearized Perturbation Dynamics')
Plot_Dynamics(t,state_lin,Full_Dynamics_Labels,'Linearized Full Dynamics')
Plot_Outputs(t,L_outputs,'Linearized Model Outputs')

% Noisy Measurements
Plot_Outputs(t,noisy_ouputs,'Noisy Measurement Model Outputs')

% Noisy States
figure();
subplot(411)
set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
plot(t,x_pert_noisy(:,1),'LineWidth',2)
hold on; grid on; grid minor
plot(t,x_pert(:,1),'--k','LineWidth',2)
xlabel('Time [s]')
ylabel(Full_Dynamics_Labels(1),'Interpreter','latex')
subplot(412)
set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
plot(t,x_pert_noisy(:,2),'LineWidth',2)
hold on; grid on; grid minor
plot(t,x_pert(:,2),'--k','LineWidth',2)
xlabel('Time [s]')
ylabel(Full_Dynamics_Labels(2),'Interpreter','latex')
subplot(413)
set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
plot(t,x_pert_noisy(:,3),'LineWidth',2)
hold on; grid on; grid minor
plot(t,x_pert(:,3),'--k','LineWidth',2)
xlabel('Time [s]')
ylabel(Full_Dynamics_Labels(3),'Interpreter','latex')
subplot(414)
set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
plot(t,x_pert_noisy(:,4),'LineWidth',2)
hold on; grid on; grid minor
plot(t,x_pert(:,4),'--k','LineWidth',2)
xlabel('Time [s]')
ylabel(Full_Dynamics_Labels(4),'Interpreter','latex')
legend('Noisy State','True State')
sgtitle('Noisy vs True States')

state_errors = x_pert_noisy - x_pert;
figure();
subplot(411)
set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
grid on; grid minor
plot(t,state_errors(:,1),'LineWidth',2)
xlabel('Time [s]')
ylabel(Error_Dynamics_Labels(1),'Interpreter','latex')
subplot(412)
set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
grid on; grid minor
plot(t,state_errors(:,2),'LineWidth',2)
xlabel('Time [s]')
ylabel(Error_Dynamics_Labels(2),'Interpreter','latex')
subplot(413)
set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
grid on; grid minor
plot(t,state_errors(:,3),'LineWidth',2)
xlabel('Time [s]')
ylabel(Error_Dynamics_Labels(3),'Interpreter','latex')
subplot(414)
set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
grid on; grid minor
plot(t,state_errors(:,4),'LineWidth',2)
xlabel('Time [s]')
ylabel(Error_Dynamics_Labels(4),'Interpreter','latex')
legend('Noisy Error')
sgtitle('Process Noise State Error')

% LKF Results
Plot_Outputs(t,LKF_outputs,'LKF Outputs')

Plot_KFState_Results(t,x_pert,x_LKF,LKF_state_err,P_LKF,'LKF State Estimate Results',...
    'LKF State Estimate Error',Full_Dynamics_Labels,Error_Dynamics_Labels)

% EKF Results
Plot_Outputs(t,EKF_outputs,'EKF Outputs')

Plot_KFState_Results(t,x_pert,x_EKF',EKF_state_err,P_EKF,'EKF State Estimate Results',...
    'EKF State Estimate Error',Full_Dynamics_Labels,Error_Dynamics_Labels)

% UKF Results
Plot_Outputs(t,UKF_outputs,'UKF Outputs')

Plot_KFState_Results(t,x_pert,x_UKF',UKF_state_err,P_UKF,'UKF State Estimate Results',...
    'UKF State Estimate Error',Full_Dynamics_Labels,Error_Dynamics_Labels)

%% Question 6

%LKF with nominal measurements
[x_LKF_log, y_LKF_log, P_LKF_log, innov_LKF_log, delta_x_LKF_log, Sv_LKF_log] = LKF(Fk, G, Hk, Q_LKF, R, Omegabar, delta_x0, Pp0, x_nom, u_nom, u, y_nom_log, ydatalog_mod);
LKF_outputs_log = [y_LKF_log, station_vis_datalog];

%LKF Plotting
% figure;
% for i = 1:4
%     subplot(4,1,i); 
%     hold on; grid on; grid minor;
% 
%     xhat_i = x_LKF_log(:,i);
%     sigma_i = sqrt(squeeze(P_LKF_log(i,i,:)));% 1σ
% 
%     plot(t_log, xhat_i, 'LineWidth', 1.5);
%     plot(t_log, xhat_i + 2*sigma_i, '--', 'LineWidth', 1); % +2σ
%     plot(t_log, xhat_i - 2*sigma_i, '--', 'LineWidth', 1); % -2σ
% 
%     ylabel(Full_Dynamics_Labels{i}, 'Interpreter','latex');
% end
% xlabel('Time [s]');
% sgtitle('LKF State Estimates with 2\sigma Bounds');

t_log = t;                  % time vector
xhat  = x_LKF_log;          % state estimates
Plog  = P_LKF_log;          % covariance logs
figure;

for i = 1:4
    
    subplot(4,2,2*i-1);
    hold on; 
    grid on; 
    grid minor;

    xhat_i  = x_LKF_log(:,i);
    sigma_i = sqrt(squeeze(P_LKF_log(i,i,:)));

    % Plot estimate + ±2σ
    p1 = plot(t, xhat_i, 'b', 'LineWidth', 1.4);
    p2 = plot(t, xhat_i + 2*sigma_i, 'r--', 'LineWidth', 1);
    p3 = plot(t, xhat_i - 2*sigma_i, 'r--', 'LineWidth', 1);

    ylabel(Full_Dynamics_Labels{i}, 'Interpreter','latex');

    % Add legend (only once per column)
    if i == 1
        legend([p1 p2], {'Estimate', '±2σ'}, 'Location','best');
    end

    title(sprintf('State %d: 2σ Envelope', i));

    subplot(4,2,2*i);
    hold on; 
    grid on; 
    grid minor;

    % Plot estimate + ±100σ
    p4 = plot(t, xhat_i, 'b', 'LineWidth', 1.4);
    p5 = plot(t, xhat_i + 100*sigma_i, 'm--', 'LineWidth', 1.2);
    p6 = plot(t, xhat_i - 100*sigma_i, 'm--', 'LineWidth', 1.2);

    ylabel(Full_Dynamics_Labels{i}, 'Interpreter','latex');

    if i == 1
        legend([p4 p5], {'Estimate', '±100σ'}, 'Location','best');
    end

    title(sprintf('State %d: 100σ (scaled) Envelope', i));
end

xlabel('Time [s]');
sgtitle('LKF State Estimates: 2σ and Scaled 100σ Envelopes');
