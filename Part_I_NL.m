clc; clear; close all

%% Simulate Non-Linear ODEs and Measured Outputs

% Constants
mu = 398600; % km^3/s^2
r0 = 6678; % km
Ydot0 = r0 * sqrt(mu/r0^3); % km/s
RE = 6378; % km
wE = (2*pi)/86400; % rad/s

% ICs
var0 = [r0 0 0 Ydot0]';
perturb_x0 = [0 0.075 0 -0.021]';
var = var0 + perturb_x0;

% Simulation time
tspan = 0:10:14000; % s

% ode45 call
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode45(@(tspan,var) OrbitEOM(tspan,var,mu),tspan,var,options);

% Get outputs
[output_var,station_vis] = MeasuredOutput(t,state,RE,wE);

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

figure();
colororder('gem12')
subplot(411)
plot(t,output_var(1:12,:),'o','LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('$\rho_{i}(t)$ [km]','Interpreter','latex')
subplot(412)
plot(t,output_var(13:24,:),'o','LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('$\dot{\rho_{i}}(t)$ [km]','Interpreter','latex')
subplot(413)
plot(t,output_var(25:36,:),'o','LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('$\phi_{i}(t)$ [km]','Interpreter','latex')
subplot(414)
plot(t,station_vis,'o','LineWidth',2)
grid on; grid minor
xlabel('Time [s]')
ylabel('Station ID Number')


