function mDotF = fuelMassFlow(UInf, altitude, aircraft, fluid)
    % function: fuelMassFlow
    % --------------------------------------------------------------------
    %   Inputs:
    %       UInf = vector of freestream velocities. [m/s]
    %       altitude = the altitude. [m].
    %       engine = structure containing the engine parameters
    %       fluid = structure containing the fluid parameters
    %   Outputs:
    %       mDotF = fuel mass flow to the engine
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TODO: write a function that calculates the required thrust of the
    %real turbofan.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function Thrust = combustorThrust(phi,UInf, altitude, engine, fluid)
    % function: combustorThrust
    % --------------------------------------------------------------------
    % This function computes the thrust of the real turbofan
    %   Inputs:
    %       UInf = vector of freestream velocities. [m/s]
    %       altitude = the altitude. [m].
    %       engine = structure containing the engine parameters
    %       fluid = structure containing the fluid parameters
    %   Outputs:
    %       Thrust = vector of the thrust produced by the engine. [N]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TODO: write a function that calculates the required thrust of the
    %real turbofan.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

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
    % get density
    [~ , ~ , ~ , rho] = atmosisa(altitude);
    % compute thrust
    T = 0.5*rho*UInf.^2*airframe.S*airframe.CD0 + ...
        2.0*airframe.W^2./(pi*airframe.e*rho*UInf.^2*airframe.b^2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

