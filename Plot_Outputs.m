function Plot_Outputs(t,outputs_cell,subplot_title)
    
    % Extract outputs and station visibility cell arrays
    outputs_cell_array = outputs_cell(:,1);
    station_vis_cell_array = outputs_cell(:,2);

    % Preallocate
    N = length(station_vis_cell_array);
    station_vis = NaN(12,N);
    outputs = NaN(36,N);

    % Create plotting arrays
    for k = 1:N
        current_station = station_vis_cell_array(k); 
        current_outputs = outputs_cell_array(k);
        current_station_vis = cell2mat(current_station);
        current_outputs_vis = cell2mat(current_outputs);
        K = height(current_station_vis);
        if isempty(current_station_vis)
            continue
        else
            for j = 1:K
                station_vis(current_station_vis(j),k) = current_station_vis(j);
                outputs(3*current_station_vis(j)-3+(1:3),k) =...
                    current_outputs_vis(3*j-3+(1:3));
            end
        end

    end  

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

