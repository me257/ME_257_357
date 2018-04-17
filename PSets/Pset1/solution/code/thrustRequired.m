function T = thrustRequired(UInf,altitude,airframe)
    % function: problemSet1_problem1a
    % --------------------------------------------------------------------
    % This function computes the required thrust.
    %   Inputs:
    %       altitude = The altitude. [m].
    %       UInf = vector of the freestream velocities. [m/s]
    %       airframe =  structure containing the necessary parameters of 
    %                   the airframe.
    %   Outputs:
    %       T = vector of the required thrust. [N]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TODO: write a function that calculates the required thrust as a
    %function of the altitude and the speed, UInf. This should use the
    %parameters specific to the HondaJet.
    
    %get density
    [~ , ~ , ~ , rho] = atmosisa(altitude);
    %compute thrust
    T = 0.5*rho*UInf.^2*airframe.S*airframe.CD0 + ...
        2.0*airframe.W^2./(pi*airframe.e*rho*UInf.^2*airframe.b^2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end

