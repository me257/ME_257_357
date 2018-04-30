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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TODO:  Calculate the temperature (T4) and outlet mass flow rates from 
    %       the combustor (mDot4). Note that any changes to the cantera gas
    %       object fluid.gas remain even after function calls. The
    %       combustion should occur via an HP reactor. That is use
    %           equilibrate(fluid.gas,'HP')
    %       where fluid.gas is the cantera phase object at T3,p3,X3. 
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
    %       The provided code sets the properties of the inlet gas. Note
    %       that the fluid object is assumed to have additional fields:
    %       fluid.gas = cantera gas object
    %       fluid.fuel.name = fuel name in cti file (e.g., 'nc12h26')
    %       fluid.fuel.nC = #number of hydrogen atoms in fuel (e.g., 12)
    %       fluid.fuel.nH = #number of hydrogen atoms in fuel (e.g., 26)
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