function Plot_Dynamics(t,state,y_labels,subplot_title)
    figure();
    subplot(411)
    plot(t,state(:,1),'LineWidth',2)
    grid on; grid minor
    xlabel('Time [s]')
    ylabel(y_labels(1),'Interpreter','latex')
    subplot(412)
    plot(t,state(:,2),'LineWidth',2)
    grid on; grid minor
    xlabel('Time [s]')
    ylabel(y_labels(2),'Interpreter','latex')
    subplot(413)
    plot(t,state(:,3),'LineWidth',2)
    grid on; grid minor
    xlabel('Time [s]')
    ylabel(y_labels(3),'Interpreter','latex')
    subplot(414)
    plot(t,state(:,4),'LineWidth',2)
    grid on; grid minor
    xlabel('Time [s]')
    ylabel(y_labels(4),'Interpreter','latex')
    sgtitle(subplot_title)
end

