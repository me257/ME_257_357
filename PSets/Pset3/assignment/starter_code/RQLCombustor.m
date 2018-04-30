function [T4,mDot4] = RQLCombustor(phi,UInf, altitude, engine, fluid)
    % function: RQLcombustor
    % --------------------------------------------------------------------
    % This function returns the output mass flow vector as a function of
    % the rerouting ratio.
    %   Inputs:
    %       phi = the overall equivalence ratio for the combustor
    %       UInf = free stream velocity [m/s]
    %       engine = structure containing the engine parameters
    %       fluid = structure containing the fluid properties including the
    %           cantera solution object
    %   Outputs:
    %       mDot4 = composition of the turbine inlet gas (in mass flow) [kg/s]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TODO:  Calculate the temperature (T4) and outlet mass flow rates from 
    %       the combustor (mDot4). Note that any changes to the cantera gas
    %       object fluid.gas remain even after function calls. The
    %       combustion should occur in three parts:
    %       1)  a call to equilibrate(fluid.gas,'HP') to simulate adiabatic
    %           combustion at T3,p3,X3 (make sure to account for the 
    %           rerouted gas in X3!)
    %       2)  Quenching and dilution. NO should be removed from the
    %           equilibrated gas, and the gas should be diluted with the
    %           rerouted gas at state 3 based on a weighted average using
    %           the rerouteRatio (beta_c).
    %       3)  The diluted gas should be equilibrated again using 
    %           equilibrate(fluid.gas,'HP'). This will be the output state
    %           of the gas.
    %       Make sure you pass in the reroute ratio. This can either be
    %       done by modifying the number of function inputs or to modify
    %       the engine object to hold this parameter.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TODO:  Calculate the inlet conditions to the combustor:
    %       T3,p3
    %       This function should be equivalent to that used for the
    %       combustor. You may find it to be useful to have this as a
    %       stand-alone function, which is shared between the two.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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