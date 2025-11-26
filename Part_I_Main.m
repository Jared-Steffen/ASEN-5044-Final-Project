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
var0 = [r0 0 0 Ydot0]';
perturb_x0 = [0 0.075 0 -0.021]';
var = var0 + perturb_x0;

% Simulation time
delta_t = 10; %s
tspan = 0:delta_t:14000; % s

% ode45 call
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode45(@(tspan,var) OrbitEOM(tspan,var,mu),tspan,var,options);
[~,state_nom] = ode45(@(tspan,var0) OrbitEOM(tspan,var0,mu),tspan,var0,options);

% Get outputs
[output_var,station_vis] = MeasuredOutput(t,state,RE,wE,true);
[output_var_nom,~] = MeasuredOutput(t,state_nom,RE,wE,false);

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
N = length(tspan);
delta_xk = zeros(4, N);
delta_xk(:,1) = perturb_x0;

for k = 1:N
    t_k = tspan(k);

    % Build H at time t_k
    for j = 1:12
        H(3*j-3+(1:3),:) = Cbar(t_k,j);
        M(3*j-3+(1:3), :) = Dbar;
        % thetaik(j,k) = thetai(t_k,j); 
    end

    % Linearized measurement at time t_k
    delta_yk(:,k) = H*delta_xk(:,k);
    yk(:,k) = output_var_nom(:,k) + delta_yk(:,k);

    % Propagate perturbation to next step (except at the last time)
    if k < N
        F = eye(4) + delta_t.*Abar(t_k);
        G = delta_t.*Bbar;                  
        delta_xk(:,k+1) = F*delta_xk(:,k);
    end

    % Determine visibility off of full nonlinear visibility
    not_vis = isnan(output_var(:,k));
    yk(not_vis, k) = NaN;
end


% Reconstruct the linearized full state
state_lin = state_nom + delta_xk.';

%% Plots
figure();
subplot(411)
plot(t,state(:,1),'LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('X Position [km]')
subplot(412)
plot(t,state(:,2),'LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('Y Position [km]')
subplot(413)
plot(t,state(:,3),'LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('X Velocity [km/s]')
subplot(414)
plot(t,state(:,4),'LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('Y Velocity [km/s]')
sgtitle('Nonlinear Dynamics')

figure();
colororder('gem12')
subplot(411)
plot(t,output_var(1:3:end,:),'o','LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('$\rho_{i}(t)$ [km]','Interpreter','latex')
subplot(412)
plot(t,output_var(2:3:end,:),'o','LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('$\dot{\rho_{i}}(t)$ [km/s]','Interpreter','latex')
subplot(413)
plot(t,output_var(3:3:end,:),'o','LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('$\phi_{i}(t)$ [radians]','Interpreter','latex')
subplot(414)
plot(t,station_vis,'o','LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('Station ID Number')
sgtitle('Nonlinear Model Ouputs')

figure();
subplot(411)
plot(tspan,delta_xk(1,:),'LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('$\delta Y$ [km]','Interpreter','latex')
subplot(412)
plot(tspan,delta_xk(2,:),'LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('$\delta X$ [km]','Interpreter','latex')
subplot(413)
plot(tspan,delta_xk(3,:),'LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('$\delta\dot{X}$ [km/s]','Interpreter','latex')
subplot(414)
plot(tspan,delta_xk(4,:),'LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('$\delta\dot{Y}$ [km/s]','Interpreter','latex')
sgtitle('Linearized Pertubation Dynamics')

figure();
subplot(411)
plot(t,state_lin(:,1),'LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('X Position [km]')
subplot(412)
plot(t,state_lin(:,2),'LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('Y Position [km]')
subplot(413)
plot(t,state_lin(:,3),'LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('X Velocity [km/s]')
subplot(414)
plot(t,state_lin(:,4),'LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('Y Velocity [km/s]')
sgtitle('Linearized Full Dynamics')

figure();
colororder('gem12')
subplot(411)
plot(t,yk(1:3:end,:),'o','LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('$\rho_{i}(t)$ [km]','Interpreter','latex')
subplot(412)
plot(t,yk(2:3:end,:),'o','LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('$\dot{\rho_{i}}(t)$ [km/s]','Interpreter','latex')
subplot(413)
plot(t,yk(3:3:end,:),'o','LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('$\phi_{i}(t)$ [radians]','Interpreter','latex')
subplot(414)
plot(t,station_vis,'o','LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('Station ID Number')
sgtitle('Linearized Model Ouputs')

% figure();
% plot(t(1:end-1),wrapToPi(yk(30,:)-thetaik(:,9)'),'o')
% hold on
% yline(pi/2)
% yline(-pi/2)