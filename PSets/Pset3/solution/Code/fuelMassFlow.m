function mDotF = fuelMassFlow(UInf, altitude, aircraft, fluid)
    % function: fuelMassFlow
    % --------------------------------------------------------------------
    % This function computes the thrust of the real turbofan
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
    %find the core mass flow
    [~ , ~ , ~ , rho0] = atmosisa(altitude);
    [~ , ~ , ~ , rhoSeaLevel] = atmosisa(0.0);
    mDot = aircraft.engine.mDotSeaLevel*rho0/rhoSeaLevel; %mass flow rate of the core gas
    Drag = thrustRequired(UInf,altitude,aircraft.airframe)/aircraft.engine.nEngines;
    
    %compute the required equivalence ratio
    phi = fzero(@(phi) Drag  ...
                - combustorThrust(phi,UInf, altitude, aircraft.engine, fluid),0.5);
    
    %compute the stoichiometric mass fraction
    nSp = nSpecies(fluid.gas);
    X3 = zeros(nSp,1);
    iFuel = speciesIndex(fluid.gas,fluid.fuel.name);
    X3(iFuel)=1.0; %stoichiometric
    iO2 = speciesIndex(fluid.gas,'o2');
    iN2 = speciesIndex(fluid.gas,'n2');
    X3(iO2)= fluid.fuel.nC+fluid.fuel.nH/4.0;
    X3(iN2)=3.76*X3(iO2);
    set(fluid.gas,'X',X3)
    Y3 = massFractions(fluid.gas);
    fSt = Y3(iFuel)/(Y3(iO2)+Y3(iN2));
    
    %compute the fuel mass flow rate (overall for two engines)
    mDotF = aircraft.engine.nEngines*phi*fSt*mDot;
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
    %Freestream conditions
    [T0 , a0 , ~ , rho0] = atmosisa(altitude);
    [~ , ~ , ~ , rhoSeaLevel] = atmosisa(0.0);
    mDot = engine.mDotSeaLevel*rho0/rhoSeaLevel;
    gamma = fluid.gamma;
    R = fluid.R;
    U0 = UInf;
    
    %get the bypass ratios
    beta = engine.bypassRatio;
    
    %make the efficiencies less verbose
    n.d=engine.efficiencies.diffuser;
    n.f=engine.efficiencies.fan;
    n.c=engine.efficiencies.compressor;
    n.b=engine.efficiencies.combustor;
    n.t=engine.efficiencies.turbine;
    n.n=engine.efficiencies.coreNozzle;
    n.n1=engine.efficiencies.fanNozzle;
    
    %get compressor ratios
    pi.f = engine.fanCompressionRatio;
    pi.c = engine.overallCompressionRatio./pi.f;
    
    %get turbine inlet temperature ratio
    [T4,~] = combustor(phi,UInf, altitude, engine, fluid); 
    tau.l = T4./T0;
    
    %compute the reference conditions
    M0 = U0./sqrt(gamma.*R.*T0);
    tau.r = 1+(gamma-1)./2.*M0.^2;
    pi.r = tau.r.^(gamma/(gamma-1.0));
    
    %diffuser pressure ratio
    tau.d = 1.0; %by assumption
    pi.d = (1.0+n.d.*(tau.r-1.0)).^(gamma./(gamma-1.0))./pi.r;
    
    %compute the fan and compressor temperature ratois
    tau.c = 1.0+1.0./n.c.*(pi.c.^((gamma-1.0)./gamma)-1.0);
    tau.f = 1.0+1.0./n.f.*(pi.f.^((gamma-1.0)./gamma)-1.0);
    
    %compute the turbine temperature ratio and pressure ratio
    tau.hpt = 1.0+tau.r.*tau.d.*tau.f.*(1.0-tau.c)./tau.l;
    tau.lpt = 1.0+tau.r.*tau.d.*(1-tau.f).*(1.0+beta)./(tau.l.*tau.hpt); %lowpressure turbine
    tau.t = tau.hpt.*tau.lpt;
    pi.t = (1.0+(n.t.*(tau.t-1.0))).^(gamma./(gamma-1));
    
    %calculate exit velocities
    Ue = a0.*sqrt(2.*n.n.*tau.l.*tau.t./(gamma-1.0).*(1.0-(pi.r.*pi.d.*pi.f.*pi.c.*pi.t).^((1-gamma)./gamma)));
    Ue1 = a0.*sqrt(2.0.*n.n1.*tau.r.*tau.f./(gamma-1.0).*(1.0-(pi.r.*pi.d.*pi.f).^((1.0-gamma)./gamma)));
    
    %calculate thrust
    Thrust = mDot.*(Ue-U0+beta.*(Ue1-U0));
    Thrust = Thrust.*engine.nEngines;
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
    %get density
    [~ , ~ , ~ , rho] = atmosisa(altitude);
    %compute thrust
    T = 0.5*rho*UInf.^2*airframe.S*airframe.CD0 + ...
        2.0*airframe.W^2./(pi*airframe.e*rho*UInf.^2*airframe.b^2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

