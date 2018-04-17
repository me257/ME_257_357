function Thrust = realTurbofanThrust(UInf, altitude, engine, fluid)
    % function: realTurbofanThrust
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
    tau.l = engine.turbineInletTemperature./T0;
    
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