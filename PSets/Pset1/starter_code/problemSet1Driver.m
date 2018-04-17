% Script: problemSet1Driver
% -------------------------------------------------------------------------
%   This script is used as the driver for problem set one. Feel free to
%   tinker with the parameters to get a feel for their sensitivities, but
%   the final submission must use the parameters as is.

clear; close all; clc

%define a structure containing the HondaJet's parameters
HondaJet.airframe.S = 17.3;                     %[m^2] , reference area
HondaJet.airframe.b = 12.12;                    %[m] , wing span
HondaJet.airframe.W = 41000;                    %[N] , weight (which is max takeoff weight)
HondaJet.airframe.CD0 = 0.01;                   %[] , zero-lift drag coefficient (or parasite drag coeff.)
HondaJet.airframe.e = 0.8;                      %[] , Oswald Efficiency
HondaJet.engine.overallCompressionRatio = 24;   %[]
HondaJet.engine.turbineInletTemperature=1600;   %[K]
HondaJet.engine.mDotSeaLevel = 8.0;             %[kg/s]
HondaJet.engine.nEngines = 2;                   %[]
HondaJet.engine.bypassRatio = 2.9;              %[]
HondaJet.engine.fanCompressionRatio = 2.0;      %[]
HondaJet.engine.efficiencies.diffuser = 0.95;   %[] 
HondaJet.engine.efficiencies.fan = 0.92;        %[]
HondaJet.engine.efficiencies.compressor = 0.87; %[]
HondaJet.engine.efficiencies.combustor = 1.00;  %[]
HondaJet.engine.efficiencies.turbine = 0.91;    %[]
HondaJet.engine.efficiencies.coreNozzle = 0.98; %[]
HondaJet.engine.efficiencies.fanNozzle = 0.98;  %[]


%define a structure containing the fluid
fluid.gamma = 1.4;
fluid.R = 287;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TODO: input a vector of altitudes and speeds
%altitudes = [x1,x2,x3]; [m]
%UInf = linspace(x1,x2)'; [m/s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create the plots for problem 1a
nAltitudes = length(altitudes);
nUInf = length(UInf);
Ts = zeros(nUInf,nAltitudes);
for iAltitude = 1:nAltitudes
    Ts(:,iAltitude)=thrustRequired(UInf,altitudes(iAltitude),HondaJet.airframe);
end
plotThrustRequired(UInf,Ts,altitudes);

%create the plots for problem 2c
nEngineModels=2;
Ts = zeros(nUInf,nEngineModels,nAltitudes);
for iAltitude = 1:nAltitudes
    Ts(:,1,iAltitude)=idealTurbojetThrust(UInf,altitudes(iAltitude),HondaJet.engine,fluid);
    Ts(:,2,iAltitude)=realTurbofanThrust(UInf,altitudes(iAltitude),HondaJet.engine,fluid);
end
plotThrusts(UInf,Ts,altitudes);


