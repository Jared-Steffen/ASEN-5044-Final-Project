function [xpkp1,Ppkp1,ymkp1] = UKF(t,mu,RE,wE,xp0,P0,Q,R,Omega,alpha,beta,kappa,y_data,station_vis)
    % Inputs:
    % t: time vector
    % mu: EArth's gravitational parameter
    % RE: radius of Earth
    % wE: rotational rate of Earth
    % xp0: initial state condition
    % P0: inital covariance matrix
    % Q: process noise covariance matrix
    % R: measurement noise covariance matrix
    % Omega: process noise mapping matrix
    % alpha: tuning parameter between [1e-4,1]
    % beta: tuning parameter(?), typically = 2
    % kappa: tuning paramater(?), typically = 0
    % y_data: noisy measurements (cell array)
    % station_vis: station visibilities

    % Outputs:
    % xpkp1: state estimates
    % Ppkp1: state covariance
    % ymkp1: predicted measurements

    % ode45 options
    options = odeset('RelTol',1e-6,'AbsTol',1e-9);

    % Initialize/Preallocate
    N = length(t);
    xpkp1 = zeros(length(P0),N);
    xpkp1(:,1) = xp0;
    n = length(xp0);
    Ppkp1 = zeros(length(P0),length(P0),N);
    Ppkp1(:,:,1) = P0;
    chik = zeros(n,2*n+1);
    chikp1 = zeros(n,2*n+1);
    wm = zeros(1,2*n+1);
    wc = zeros(1,2*n+1);
    ymkp1 = cell(N,1);

    % Calulate necessary constants
    lambda = alpha^2*(n+kappa)-n;
    wm(1:2*n) = 1/(2*(n+lambda));
    wm(2*n+1) = lambda/(n+lambda);
    wc(1:2*n) = wm(1:2*n);
    wc(2*n+1) = wm(2*n+1)+1-alpha^2+beta;

    % UKF Algorithm
    for k = 1:N-1

        % Step 1: Dynamics Prediction 
        % Generate 2n+1 Sigma Points for Dynamics Prediciton
        chik(:,2*n+1) = xpkp1(:,k); % mean in 9th position
        Sk = chol(Ppkp1(:,:,k),'lower');
        for i = 1:n
            chik(:,i) = xpkp1(:,k) + (sqrt(n+lambda))*Sk(:,i); % 1-4 of 2nd moment
        end
        for i = n+1:2*n
            chik(:,i) = xpkp1(:,k) - (sqrt(n+lambda))*Sk(:,i-4); % 5-8 of 2nd moment
        end
        
        % Propogate Sigma Points through NL Dynamics, get prior info
        current_tspan = [t(k) t(k+1)];
        for i = 1:2*n+1
            current_chi = chik(:,i);
            [~,chik_prop] = ode45(@(current_tspan,current_chi)...
                OrbitEOM(current_tspan,current_chi,mu),current_tspan,current_chi,options);
            chibarkp1(:,i) = chik_prop(end,:)'; 
        end
        xmkp1 = chibarkp1 * wm';
        Pmkp1 = zeros(size(P0));
        for i = 1:2*n+1
            Pmkp1 = wc(i)*(chibarkp1(:,i)-xmkp1)*(chibarkp1(:,i)-xmkp1)' + Pmkp1;
        end
        Pmkp1 = Pmkp1 + Omega*Q*Omega';

        % Step 2: Measurement Update
        % Generate 2n+1 Sigma Points for Measurement Update
        chikp1(:,2*n+1) = xmkp1; % pred mean in 9th position
        Sbarkp1 = chol(Pmkp1,'lower');
        for i = 1:n
            chikp1(:,i) = xmkp1 + (sqrt(n+lambda))*Sbarkp1(:,i); % 1-4 of 2nd moment
        end
        for i = n+1:2*n
            chikp1(:,i) = xmkp1 - (sqrt(n+lambda))*Sbarkp1(:,i-4); % 5-8 of 2nd moment
        end
        
        % Propogate Sigma Points through NL Measurements, get prior info
        for i = 1:2*n+1
            current_chikp1 = chikp1(:,i);
            [output_var,~] = MeasuredOutput...
                (t(k+1),current_chikp1',RE,wE,false);
            gammabar_unfiltered(:,i) = output_var;
        end
        vis_stations = station_vis{k+1};

        % Filter Visible Stations
        gammabar_kp1 = [];
        for j = 1:length(vis_stations)
            gammabar_kp1(3*j-3+(1:3),:) = gammabar_unfiltered(3*vis_stations(j)-3+(1:3),:);
        end
        if isempty(gammabar_kp1)
            ymkp1{k+1} = [];
            xpkp1(:,k+1) = xmkp1;
            Ppkp1(:,:,k+1) = Pmkp1;
        else
            ymkp1{k+1} = gammabar_kp1 * wm';
            K = length(ymkp1{k+1});
            Pyykp1 = zeros(size(kron(eye(K/3),R)));
            for i = 1:2*n+1
                Pyykp1 = wc(i)*(gammabar_kp1(:,i)-ymkp1{k+1})...
                    *(gammabar_kp1(:,i)-ymkp1{k+1})' + Pyykp1;
            end
            Pyykp1 = Pyykp1 + kron(eye(K/3),R);
    
            % Get state-measurement cross covariance matrix
            Pxykp1 = zeros(n,size(kron(eye(K/3),R),1));
            for i = 1:2*n+1
                Pxykp1 = wc(i)*(chibarkp1(:,i)-xmkp1)...
                    *(gammabar_kp1(:,i)-ymkp1{k+1})' + Pxykp1;
            end
    
            % Estimate Kalman Gain
            Kkp1 = Pxykp1/Pyykp1;
    
            % Update State and Covariace
            xpkp1(:,k+1) = xmkp1+Kkp1*(y_data{k+1}-ymkp1{k+1});
            Ppkp1(:,:,k+1) = Pmkp1-Kkp1*Pyykp1*Kkp1';
        end
    end

end

