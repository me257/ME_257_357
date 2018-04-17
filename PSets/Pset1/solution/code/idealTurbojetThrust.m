function Thrust = idealTurbojetThrust(UInf, altitude, engine, fluid)
    % function: idealTurbojetThrust
    % --------------------------------------------------------------------
    % This function computes the thrust of the ideal turbojet
    %   Inputs:
    %       UInf = vector of freestream velocities. [m/s]
    %       altitude = the altitude. [m].
    %       engine = structure containing the engine parameters
    %       fluid = structure containing the fluid parameters
    %   Outputs:
    %       Thrust = vector of the thrust produced by the engine. [N]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TODO: write a function that calculates the required thrust of the
    %ideal turbojet.
    
    %Freestream conditions
    [T_0 , ~ , P_0 , rho_0] = atmosisa(altitude);
    [~ , ~ , ~ , rhoSeaLevel] = atmosisa(0.0);
    mDot = engine.mDotSeaLevel*rho_0/rhoSeaLevel;
    cp = fluid.R*fluid.gamma/(fluid.gamma-1.0);

    nUInf = length(UInf);
    T = zeros(nUInf,9);
    P = zeros(nUInf,9);
    U = zeros(nUInf,9);

    %Freestream to inlet opening: 0 -> 1
    T(:,1) = T_0;
    P(:,1) = P_0;
    U(:,1) = UInf;

    %Inlet diffuser: 1 -> 2
    T(:,2) = T(:,1) + 1/(2*cp) * U(:,1).^2;
    P(:,2) = P(:,1) .* ( T(:,2) ./ T(:,1) ) .^ (fluid.gamma/(fluid.gamma-1));


    %Compressor: 2 -> 3
    P(:,3) = P(:,2) * engine.overallCompressionRatio;
    T(:,3) = T(:,2) .* ( P(:,3) ./ P(:,2) ) .^ ((fluid.gamma-1)/fluid.gamma);


    %Combustor: 3 -> 4
    P(:,4) = P(:,3);
    T(:,4) = engine.turbineInletTemperature;


    %Turbine: 4 -> 5
    T(:,5) = T(:,4) - ( T(:,3) - T(:,2) );
    P(:,5) = P(:,4) .* ( T(:,5) ./ T(:,4) ) .^ (fluid.gamma/(fluid.gamma-1));


    %Afterburner: 5 -> 6
    T(:,6) = T(:,5);
    P(:,6) = P(:,5);

    %No losses 6 -> 7
    T(:,7) = T(:,6);
    P(:,7) = P(:,6);

    %Exhaust nozzle throat: 7 -> 8
    T(:,8) = T(:,7) ./ ( 1 + (fluid.gamma-1)/2 );
    P(:,8) = P(:,7) .* ( T(:,8) ./ T(:,7) ) .^ (fluid.gamma/(fluid.gamma-1));
    U(:,8) = sqrt(fluid.gamma * fluid.R * T(:,8));

    %Nozzle throat to exit: 8 -> e or 9
    P(:,9) = P_0;
    T(:,9) = T(:,8) .* ( P(:,9) ./ P(:,8) ) .^ ((fluid.gamma-1)/fluid.gamma);
    U(:,9) = sqrt( 2*cp*(T(:,8) - T(:,9)) + U(:,8).^2);

    Thrust = mDot * ( U(:,9) - UInf );
    Thrust = Thrust*engine.nEngines;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end