function [output_var,station_vis] = MeasuredOutput(t,state,RE,wE)
    % Goal: Output measured outputs (rho, rho_dot, phi) for each visible station

    % Extract state info
    x_t = state(:,1); % km
    x_dot_t = state(:,2); % km/s
    y_t = state(:,3); % km
    y_dot_t = state(:,4); % km/s
    
    % Get positions and viewing angle for each station
    for i = 1:12
        theta_i0 = (i-1)*pi/6; % rad
        Xs_it = RE*cos(wE.*t+theta_i0); % km
        Ys_it = RE*sin(wE.*t+theta_i0); % km
        Xsdot_it = -RE*wE*sin(wE.*t+theta_i0); % km/s
        Ysdot_it = RE*wE*cos(wE.*t+theta_i0); % km/s
        theta_it = atan2(Ys_it,Xs_it); % rad
    
        % Determine visibility and outputs
        phi_it(:,i) = atan2(y_t-Ys_it,x_t-Xs_it);
        rho_it(:,i) = sqrt((x_t-Xs_it).^2 + (y_t-Ys_it).^2); 
        rhodot_it(:,i) = (((x_t - Xs_it).*(x_dot_t - Xsdot_it)) +...
            ((y_t-Ys_it).*(y_dot_t-Ysdot_it)))./rho_it(:,i);
        station_vis(:,i) = i.*ones(length(rho_it(:,i)),1);
        for k = 1:length(t)
            if abs(wrapToPi(phi_it(k,i) - theta_it(k))) < pi/2
                phi_it(k,i) = phi_it(k,i);
                rho_it(k,i) = rho_it(k,i);
                rhodot_it(k,i) = rhodot_it(k,i);
                station_vis(k,i) = station_vis(k,i);
            else
                phi_it(k,i) = NaN;
                rho_it(k,i) = NaN;
                rhodot_it(k,i) = NaN;
                station_vis(k,i) = NaN;
            end
        end
    
    end
    
    output_var = [rho_it';rhodot_it';phi_it'];
end