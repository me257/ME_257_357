function plotThrusts(UInf,Ts,altitudes)
    % function: plotThrusts
    % --------------------------------------------------------------------
    % This function builds the plots required for problem 2.c.i of problem 
    % set 1.
    %   Inputs:
    %       UInf = vector of freestream velocities. [m/s]
    %       Ts = an array of the thrusts (nUInf X nEngines X nAltitudes).
    %       [N]
    %       altitudes=vector of altitudes [m]. 
    %   Outputs:
    %       plots of thrust vs UInf for the inputed altitudes and the two
    %       engines.
    
    %create the plots
    figure()
    nEngineModels=2;
    legendEntries = cell(1, length(altitudes)*nEngineModels);
    index = @(i,j) (i-1)*length(altitudes)+j;
    engines = {'turobojet','turbofan'};
    colors = zeros(length(altitudes),3);
    colors(:,1)=((1:length(altitudes))-1)/(length(altitudes)-1);
    linestyles = {'-','--'};
    for i = 1:nEngineModels
        for j = 1:length(altitudes)	
                plot(UInf, Ts(:,i,j),linestyles{i},'color', colors(j,:))
                hold on
                legendEntries{1,index(i,j)} = [engines{i}, ...
                                           ': altitude = ', ... 
                                           num2str(altitudes(j)),  ' m'];
        end
    end
    hold off
    title('Problem 2.c.i')
    xlabel('$U_\infty$ [m/s]')
    ylabel('Thrust [N]')
    legend(legendEntries,'location','best')

end