function Plot_Outputs(t,outputs,station_vis,subplot_title)
    figure();
    colororder('gem12')
    subplot(411)
    plot(t,outputs(1:3:end,:),'o','LineWidth',2)
    grid on; grid minor
    xlabel('Time [s]')
    ylabel('$\rho_{i}(t)$ [km]','Interpreter','latex')
    subplot(412)
    plot(t,outputs(2:3:end,:),'o','LineWidth',2)
    grid on; grid minor
    xlabel('Time [s]')
    ylabel('$\dot{\rho_{i}}(t)$ [km/s]','Interpreter','latex')
    subplot(413)
    plot(t,outputs(3:3:end,:),'o','LineWidth',2)
    grid on; grid minor
    xlabel('Time [s]')
    ylabel('$\phi_{i}(t)$ [radians]','Interpreter','latex')
    subplot(414)
    plot(t,station_vis,'o','LineWidth',2)
    grid on; grid minor
    xlabel('Time [s]')
    ylabel('Station ID Number')
    sgtitle(subplot_title)
end

