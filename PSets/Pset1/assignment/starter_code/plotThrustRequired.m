function plotThrustRequired(UInf,Ts,altitudes)
    % function: problemSet1_problem1a
    % --------------------------------------------------------------------
    % This function builds the plots required for problem 1a of problem set
    % 1.
    %   Inputs:
    %       UInf = vector of freestream velocities. [m/s]
    %       Ts = array of required thrust at different altitudes. [N]
    %       altitudes=vector of altitudes [m]. 
    %   Outputs:
    %       plots of required thrust vs UInf for the inputed altitudes.

    %This loop builds the required thrust vector at the given altitudes
    MInf = zeros(length(UInf), length(altitudes));
    nAltitudes = length(altitudes);
    for i = 1:nAltitudes
        %atmosisa provides the thermodynamic data as a function of altitude
        [~ , a , ~ , ~] = atmosisa(altitudes(i,1));
        MInf(:,i) = UInf / a;
    end
    
    %this loop creates the strings used for the legend
    legendEntries = cell(1, nAltitudes);
    for i = 1:nAltitudes	
        legendEntries{1,i} = ['altitude = ' num2str(altitudes(i)) ' m'];
    end
    
    %create the plots
    figure()
    plot(UInf, Ts)
    hold on
    title('Problem 1.a')
    xlabel('U_\infty [m/s]')
    ylabel('Required Thrust [N]')
    legend(legendEntries)

    figure()
    plot(MInf, Ts, '-')
    hold on
    title('Problem 1.a')
    xlabel('M_\infty')
    ylabel('Required Thrust [N]')
    legend(legendEntries)

end