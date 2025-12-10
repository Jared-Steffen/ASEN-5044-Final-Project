function Plot_KFState_Results(t,state,estimated_state,state_errors,Ppkp1,...
    subplot_title1,subplot_title2,y_labels1,y_labels2)

    % subplot_title1 and y_labels1 are for the state estimate
    % subplot_title2 and y_labels2 are for the estimate error
    
    % Extract Variances
    X_var     = squeeze(Ppkp1(1,1,:));
    X_dot_var = squeeze(Ppkp1(2,2,:));
    Y_var     = squeeze(Ppkp1(3,3,:));
    Y_dot_var = squeeze(Ppkp1(4,4,:));

    % Calculate Sigma Bound Magnitude
    X_sigma(:,1) = sqrt(X_var);
    X_dot_sigma(:,1) = sqrt(X_dot_var);
    Y_sigma(:,1) = sqrt(Y_var);
    Y_dot_sigma(:,1) = sqrt(Y_dot_var);

    figure();
    subplot(411)
    set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
    plot(t,estimated_state(:,1),'LineWidth',2)
    hold on; grid on; grid minor
    plot(t,state(:,1),'--k','LineWidth',2)
    xlabel('Time [s]')
    ylabel(y_labels1(1),'Interpreter','latex')
    subplot(412)
    set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
    plot(t,estimated_state(:,2),'LineWidth',2)
    hold on; grid on; grid minor
    plot(t,state(:,2),'--k','LineWidth',2)
    xlabel('Time [s]')
    ylabel(y_labels1(2),'Interpreter','latex')
    subplot(413)
    set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
    plot(t,estimated_state(:,3),'LineWidth',2)
    hold on; grid on; grid minor
    plot(t,state(:,3),'--k','LineWidth',2)
    xlabel('Time [s]')
    ylabel(y_labels1(3),'Interpreter','latex')
    subplot(414)
    set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
    plot(t,estimated_state(:,4),'LineWidth',2)
    hold on; grid on; grid minor
    plot(t,state(:,4),'--k','LineWidth',2)
    xlabel('Time [s]')
    ylabel(y_labels1(4),'Interpreter','latex')
    legend('Estimated State','True State')
    sgtitle(subplot_title1)

    figure();
    subplot(411)
    set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
    hold on; grid on; grid minor
    plot(t,state_errors(:,1),'LineWidth',2)
    plot(t,2*X_sigma,':','LineWidth',2)
    plot(t,-2*X_sigma,':','LineWidth',2)
    xlabel('Time [s]')
    ylabel(y_labels2(1),'Interpreter','latex')
    subplot(412)
    set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
    hold on; grid on; grid minor
    plot(t,state_errors(:,2),'LineWidth',2)
    plot(t,2*X_dot_sigma,':','LineWidth',2)
    plot(t,-2*X_dot_sigma,':','LineWidth',2)
    xlabel('Time [s]')
    ylabel(y_labels2(2),'Interpreter','latex')
    subplot(413)
    set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
    hold on; grid on; grid minor
    plot(t,state_errors(:,3),'LineWidth',2)
    plot(t,2*Y_sigma,':','LineWidth',2)
    plot(t,-2*Y_sigma,':','LineWidth',2)
    xlabel('Time [s]')
    ylabel(y_labels2(3),'Interpreter','latex')
    subplot(414)
    set(gca, 'ColorOrder', [0 0.4470 0.7410], 'NextPlot', 'add'); 
    hold on; grid on; grid minor
    plot(t,state_errors(:,4),'LineWidth',2)
    plot(t,2*Y_dot_sigma,':','LineWidth',2)
    plot(t,-2*Y_dot_sigma,':','LineWidth',2)
    xlabel('Time [s]')
    ylabel(y_labels2(4),'Interpreter','latex')
    legend('Estimated Error','\pm 2\sigma Bounds','')
    sgtitle(subplot_title2)
    
end

