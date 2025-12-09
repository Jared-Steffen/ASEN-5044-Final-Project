clc; clear all; close all;

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
N = length(tspan);

% Initialize Covariance
Pp0 = diag([1000,1,1000,1]);

% ode45 options
options = odeset('RelTol',1e-6,'AbsTol',1e-9);

% Loading Qtrue and Rtrue
data = load('orbitdeterm_finalproj_KFdata.mat');
Qtrue = data.Qtrue;
R = data.Rtrue;
% KF Tuning Knob
Q = 80*Qtrue;

% Inputs
u_nom = zeros(2,length(tspan));
u = zeros(2,length(tspan));

% Getting nominal trajectory and measurements
[~,x_nom] = ode45(@(tspan,var0) OrbitEOM(tspan,var0,mu),tspan,nom_var0,options);
[output_var_nom,~] = MeasuredOutput(tspan',x_nom,RE,wE,false);

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

Gamma = Bbar; % w1,2 has same mapping as u1,2
Omegabar = delta_t.*Gamma;

% Compute Fk once outside MC loop
Fk = zeros(4,4,N);
for k = 1:N-1
    Fk(:,:,k) = eye(4) + delta_t.*Abar(tspan(k));
end
G = delta_t.*Bbar;

%% MC
num_sims = 10;
epsilon_x = zeros(N,num_sims);
epsilon_y = zeros(N,num_sims);
for sim_num = 1:num_sims
    %% Simulating nonlinear EOMs with perturbed initial condition and process noise
    % ode45 call - pert w/ noise
    x_true_pert = pert_var0;
    x_pert_noisy = zeros(N,4);
    x_pert_noisy(1,:) = x_true_pert';
    e_x_EKF = zeros(N,4);
    for k = 1:N-1
        % Propagate perturbed true state with noise using ode45
        current_tspan = [tspan(k) tspan(k+1)];
        [~,x_pert_k] = ode45(@(current_tspan,pert_var0)...
            OrbitEOM(current_tspan,pert_var0,mu),current_tspan,x_true_pert,options);
        xtrue_kp1 = x_pert_k(end,:); 
        w_k = chol(Qtrue,"lower")*randn(2,1);
        x_true_pert = xtrue_kp1' + Omegabar*w_k;
        x_pert_noisy(k+1,:) = x_true_pert;
    end

    %% Simulating measurements for noisy NL simulation, adding measurment noise
    % Also, building H matrices at each k time step
    [station_meas,station_vis] = MeasuredOutput(tspan',x_pert_noisy,RE,wE,true);
    % Convert into cell array - filter nominal measurements to match NL idxs
    y_nom = cell(N,1);
    y_pert_noisy = cell(N,1);
    station_vis_cell = cell(N,1);

    delta_xk = zeros(4, N);
    delta_xk(:,1) = delta_x0;
    Hk = cell(N,1);
    p(1) = 0;
    for k = 2:N
        y_pert = station_meas(~any(isnan(station_meas(:,k)),2), k);
        y_nom{k} = output_var_nom(~any(isnan(station_meas(:,k)),2), k);
        K = length(y_pert);
        Sv = chol(kron(eye(K/3),R),'lower');
        qk = randn(K,1);
        y_pert_noisy{k} = y_pert + Sv*qk;
        station_vis_cell{k} = station_vis(~any(isnan(station_vis(:,k)),2), k);

        % Build H at time t_k
        num_vis = station_vis_cell{k};
        H = zeros(3*length(num_vis),4);
        for i = 1:length(num_vis)
            H(3*i-3+(1:3),:) = Cbar(tspan(k),num_vis(i));
        end
        Hk{k} = H;

        % p keeps track of the number of available measurements at step k
        p(k) = size(H,1);
    end
    
    % Run LKF and calculate NEES and NIS
    delta_x0_KF = zeros(4,1);
    [x_LKF,y_LKF,Pmkp1_LKF,Ppkp1_LKF,innov_LKF,delta_x_LKF,Skp1] = LKF...
        (Fk,G,Hk,Q,R,Omegabar,delta_x0_KF,Pp0,x_nom,u_nom,u,y_nom,y_pert_noisy);
    % e_y_LKF = cell(N,1);
    for k=2:N
        % calculate epsilon_x for each time step k
        e_x_LKF(k,:) = x_pert_noisy(k,:) - x_LKF(k,:);
        epsilon_x_LKF(k,sim_num) = e_x_LKF(k,:) * (Ppkp1_LKF(:,:,k) \ e_x_LKF(k,:)');
        % calculate epsilon_y for each time step k
        epsilon_y_LKF(k,sim_num) = innov_LKF{k}' * (Skp1{k} \ innov_LKF{k});
    end

    % Run EKF and calculate NEES and NIS
    [x_EKF, P_EKF, y_EKF, Sv_EKF] = EKF...
        (Q, R, y_pert_noisy, tspan', mu, RE, wE, nom_var0, Pp0, station_vis_cell);
    e_y_EKF = cell(N,1);
    for k=2:N
        % calculate epsilon_x for each time step k
        e_x_EKF(k,:) = x_pert_noisy(k,:) - x_EKF(:,k)';
        epsilon_x_EKF(k,sim_num) = e_x_EKF(k,:) * (P_EKF(:,:,k) \ e_x_EKF(k,:)');
        % calculate epsilon_y for each time step k
        e_y_EKF{k} = y_pert_noisy{k} - y_EKF{k};
        epsilon_y_EKF(k,sim_num) = e_y_EKF{k}' * (Sv_EKF{k} \ e_y_EKF{k});
    end
end
% Calculate NEES and NIS
epsilonbar_x_LKF = (1/num_sims) .* sum(epsilon_x_LKF,2);
epsilonbar_y_LKF = (1/num_sims) .* sum(epsilon_y_LKF,2);

% Getting upper and lower bounds for NEES and NIS
alpha = 0.05;
r1_NEES = chi2inv(alpha/2,num_sims*4)./num_sims;
r2_NEES = chi2inv(1-alpha/2,num_sims*4)./num_sims;
r1_NIS = chi2inv(alpha/2,num_sims*p)./num_sims;
r2_NIS = chi2inv(1-alpha/2,num_sims*p)./num_sims;

%% Plots of NEES and NIS for LKF and EKF
figure;
scatter(1:N, epsilonbar_x_LKF, 16, 'filled'); grid on; hold on;
xlabel('Time Step k');
ylabel('$\bar{\epsilon}_x$', 'Interpreter', 'latex');
title('LKF NEES Test Results');
xlim([1,N]);
yline(r1_NEES,'--',"r1","Color","red");
yline(r2_NEES,'--',"r2","Color","red");

figure;
scatter(1:N, epsilonbar_y_LKF, 16, 'filled'); grid on; hold on;
xlabel('Time Step k');
ylabel('$\bar{\epsilon}_y$', 'Interpreter', 'latex');
title('LKF NIS Test Results');
xlim([1,N]);
plot(1:N,r1_NIS,'--',"Color","red");
plot(1:N,r2_NIS,'--',"Color","red");

epsilonbar_x_EKF = (1/num_sims) .* sum(epsilon_x_EKF,2);
epsilonbar_y_EKF = (1/num_sims) .* sum(epsilon_y_EKF,2);

figure;
scatter(1:N, epsilonbar_x_EKF, 16, 'filled'); grid on; hold on;
xlabel('Time Step k');
ylabel('$\bar{\epsilon}_x$', 'Interpreter', 'latex');
title('EKF NEES Test Results');
xlim([1,N]);
yline(r1_NEES,'--',"r1","Color","red");
yline(r2_NEES,'--',"r2","Color","red");

figure;
scatter(1:N, epsilonbar_y_EKF, 16, 'filled'); grid on; hold on;
xlabel('Time Step k');
ylabel('$\bar{\epsilon}_y$', 'Interpreter', 'latex');
title('EKF NIS Test Results');
xlim([1,N]);
% plot(1:N,r1_NIS,'--',"Color","red");
% plot(1:N,r2_NIS,'--',"Color","red");