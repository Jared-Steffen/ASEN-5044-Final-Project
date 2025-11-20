function [var_dot] = OrbitEOM(~,var,mu)
    % Goal: Output ODEs for ode45 where var is a 4x1 state vector (2D pos and vel)

    % Extract state variables
    x = var(1);
    u = var(2);
    y = var(3);
    v = var(4);

    % Calculate radius
    r = sqrt(x^2 + y^2);

    % Assign t.r.o.c variables
    x_dot = u;
    y_dot = v;
    u_dot = (-mu*x)/r^3;
    v_dot = (-mu*y)/r^3;

    % Final state derivative
    var_dot = [x_dot;u_dot;y_dot;v_dot];
end