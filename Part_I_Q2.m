clc; clear all; close all;

% Constants
mu = 398600; % km^3/s^2
r0 = 6678; % km
w = sqrt(mu/r0^3);
RE = 6378; % km
wE = (2*pi)/86400; % rad/s
delta_T = 10; % s

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
Xidot = @(t,i) -wE * RE * sin(wE*t + thetai0(i));
Yidot = @(t,i) wE * RE * cos(wE*t + thetai0(i));

delta_Xnom = @(t,i) r0*cos(w*t) - Xi(t,i);
delta_Ynom = @(t,i) r0*sin(w*t) - Yi(t,i);
delta_Xdotnom = @(t,i) -w*r0*sin(w*t) - Xidot(t,i);
delta_Ydotnom = @(t,i) w*r0*cos(w*t) - Yidot(t,i);
Anom = @(t,i) delta_Xnom(t,i)*delta_Xdotnom(t,i)...
    + delta_Ynom(t,i)*delta_Ydotnom(t,i);
rhonom = @(t,i) sqrt(delta_Xnom(t,i)^2 + delta_Ynom(t,i)^2);

Cbar = @(t,i) [delta_Xnom(t,i)/rhonom(t,i), 0, delta_Ynom(t,i)/rhonom(t,i), 0 % row 1
(rhonom(t,i)*delta_Xdotnom(t,i) - Anom(t,i)*(delta_Xnom(t,i)/rhonom(t,i))) / rhonom(t,i)^2,... row 2 col 1
delta_Xnom(t,i)/rhonom(t,i),... row 2 col 2
(rhonom(t,i)*delta_Ydotnom(t,i) - Anom(t,i)*(delta_Ynom(t,i)/rhonom(t,i))) / rhonom(t,i)^2,... row 2 col 3
delta_Ynom(t,i)/rhonom(t,i) % row 2 col 4
-delta_Ynom(t,i)/rhonom(t,i)^2, 0, delta_Xnom(t,i)/rhonom(t,i)^2, 0]; % row 4

Dbar = zeros(3,2);