function [T4,mDot4] = combustor(phi,UInf, altitude, engine, fluid)
    % function: combustor
    % --------------------------------------------------------------------
    % This function returns the turbine inlet temperature and mass
    % fractions.
    %   Inputs:
    %       phi = equivalence ratio
    %       UInf = free stream velocity [m/s]
    %       engine = structure containing the engine parameters
    %       fluid = structure containing the fluid properties including the
    %           cantera solution object
    %   Outputs:
    %       T4 = The turbine inlet temperature [K]
    %       mDot4 = composition of the turbine inlet gas (in mass flow) [kg/s]
    
    %set the fluid.gas object to the inlet conditions
    setCombustorInletConditions(phi, UInf, altitude, engine, fluid);
    
    %determine the mass flow
    [~ , ~ , ~ , rho0] = atmosisa(altitude);
    [~ , ~ , ~ , rhoSeaLevel] = atmosisa(0.0);
    mDotA = engine.mDotSeaLevel*rho0/rhoSeaLevel;
    mDotC = mDotA/(1+engine.bypassRatio);
    f = fuelAirRatio(fluid);
    mDotC = (1+f)*mDotC;
    
    %combust the fuel and get mass flow of the species
    equilibrate(fluid.gas,'HP');
    T4 = temperature(fluid.gas);
    mDot4 = massFractions(fluid.gas)*mDotC;
end
function setCombustorInletConditions(phi,UInf, altitude, engine, fluid)
    % function: setCombustorInletConditions
    % --------------------------------------------------------------------
    % This function sets the fluid.gas solution object to the combustor
    % inlet conditions after the fuel is mixed with the oxidizer
    %   Inputs:
    %       phi = equivalence ratio
    %       UInf = free stream velocity [m/s]
    %       engine = structure containing the engine parameters
    %       fluid = structure containing the fluid properties including the
    %           cantera solution object
    %   Outputs:
    %       fluid.gas object is altered to inlet state
    %Freestream conditions
    [T0 , ~ , p0 , ~] = atmosisa(altitude);
    gamma = fluid.gamma;
    R = fluid.R;
    U0 = UInf; 
    
    %make the efficiencies less verbose
    n.d=engine.efficiencies.diffuser;
    n.f=engine.efficiencies.fan;
    n.c=engine.efficiencies.compressor;
    n.b=engine.efficiencies.combustor;
    
    %get compressor ratios
    pi.f = engine.fanCompressionRatio;
    pi.c = engine.overallCompressionRatio./pi.f;
    
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
    
    %get the combustor inlet conditions (assume M4<<1)
    p3 = p0*pi.r*pi.d*pi.f*pi.c;
    T3 = T0*tau.r*tau.d*tau.f*tau.c;
    
    %set the inlet gas
    nSp = nSpecies(fluid.gas);
    X3 = zeros(nSp,1);
    iFuel = speciesIndex(fluid.gas,fluid.fuel.name);
    X3(iFuel)=phi;
    iO2 = speciesIndex(fluid.gas,'o2');
    iN2 = speciesIndex(fluid.gas,'n2');
    X3(iO2)= fluid.fuel.nC+fluid.fuel.nH/4.0;
    X3(iN2)=3.76*X3(iO2);
    set(fluid.gas,'T',T3,'P',p3,'X',X3)
end

function f = fuelAirRatio(fluid)
    % function: fuelAirRatio
    % --------------------------------------------------------------------
    % Determines the fuel-air ratio of the unburned gas.
    %   Inputs:
    %       fluid = structure containing the fluid properties including the
    %           cantera solution object
    %   Outputs:
    %       f = fuel-air ratio
    Y = massFractions(fluid.gas);
    iFuel = speciesIndex(fluid.gas,fluid.fuel.name);
    iO2 = speciesIndex(fluid.gas,'o2');
    iN2 = speciesIndex(fluid.gas,'n2');
    f = Y(iFuel)/(Y(iO2)+Y(iN2));
end