% Script: problemSet1Driver
% -------------------------------------------------------------------------
%   This script is used as the driver for problem set one. It simply
%   requires running and no changes are necessary. Feel free to tinker with
%   the parameters herein to get a feel of their effects not explicity
%   asked within the problem set.

clear; close all; clc

%constants
ftToM = 12*0.0254;
g=9.81;

%define a structure containing the HondaJet's parameters
HondaJet.airframe.S = 17.3;                     %[m^2] , reference area
HondaJet.airframe.b = 12.12;                    %[m] , wing span
HondaJet.airframe.W = 41000;                    %[N] , weight (which is max takeoff weight)
HondaJet.airframe.CD0 = 0.01;                   %[] , zero-lift drag coefficient (or parasite drag coeff.)
HondaJet.airframe.e = 0.8;                      %[] , Oswald Efficiency
HondaJet.engine.overallCompressionRatio = 24;           %[]
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
HondaJet.engine.mDotNO = 0.04e-3;               %[kg/s]
HondaJet.engine.combustorRerouteRatio = 0.0;    %[]
HondaJet.airframe.WFuel = 5600;                 %[N] 
HondaJet.airframe.WInitial = HondaJet.airframe.W;%[N]


%define a structure containing the fluid
fluid.gas = Solution('ndodecane_mech.xml');
fluid.fuel.name = 'nc12h26';
fluid.fuel.nC = 12;
fluid.fuel.nH = 26;
fluid.gamma = 1.4;
fluid.R = 287;
iNO = speciesIndex(fluid.gas,'no');
hs = [0,30000,43000]*ftToM; %[m]
UInfs = [90,150,200]; %[m/s]

nPhi = 100;
phis = logspace(log10(0.15),log10(5),nPhi);
T4s = zeros(nPhi,1);
mNOs = zeros(nPhi,1);
thrusts = zeros(nPhi,1);
nH = length(hs);
labels = cell(nH,1);

%% set up figures
temperatureFigure = figure();
thrustFigure = figure();

for iH = 1:nH
    T4Original = HondaJet.engine.turbineInletTemperature;
    for iPhi = 1:nPhi
        [T4,mDot4] = combustor(phis(iPhi),UInfs(iH), hs(iH), HondaJet.engine, fluid);
        HondaJet.engine.turbineInletTemperature=T4;
        thrusts(iPhi)=realTurbofanThrust(UInfs(iH), hs(iH), HondaJet.engine, fluid);
        mNOs(iPhi)=mDot4(iNO);
        T4s(iPhi)=T4;
    end
    HondaJet.engine.turbineInletTemperature=T4Original;
    labels{iH}=['h=',num2str(hs(iH)),' m'];
    
    figure(thrustFigure);
    semilogx(phis,thrusts)
    hold on
    
    figure(temperatureFigure);
    semilogx(phis,T4s/HondaJet.engine.turbineInletTemperature)
    hold on
end

figure(thrustFigure);
xlabel('$\phi$')
ylabel('$T\ [\mathrm{N}]$')
legend(labels,'location','best')
hold off
figure(temperatureFigure);
xlabel('$\phi$')
ylabel('$T_4/T_{4,\mathrm{max}}$')
legend(labels,'location','best')
hold off


%% plot the mNOs figures
mDotNOFigure = figure();
for iPhi = 1:nPhi
        [T4,mDot4] = combustor(phis(iPhi),UInfs(1), hs(1), HondaJet.engine, fluid);
        thrusts(iPhi)=realTurbofanThrust(UInfs(1), hs(1), HondaJet.engine, fluid);
        mNOs(iPhi)=mDot4(iNO);
        T4s(iPhi)=T4;
end
loglog(phis,mNOs/HondaJet.engine.mDotNO)
hold on
betaC = [1,5,10];
nBetaC = length(betaC);
labels = cell(nBetaC+1,1);
labels{1}='$\beta_\mathrm{C}$=0';
for iBetaC = 1:nBetaC
    HondaJet.engine.combustorRerouteRatio = betaC(iBetaC);
    for iPhi = 1:nPhi
        [T4,mDot4] = RQLCombustor(phis(iPhi),UInfs(1), hs(1), HondaJet.engine, fluid);
        thrusts(iPhi)=realTurbofanThrust(UInfs(1), hs(1), HondaJet.engine, fluid);
        mNOs(iPhi)=mDot4(iNO);
        T4s(iPhi)=T4;
    end
    loglog(phis,mNOs/HondaJet.engine.mDotNO)
    labels{iBetaC+1}=['$\beta_\mathrm{C}$=',num2str(betaC(iBetaC))];
end

figure(mDotNOFigure);
xlabel('$\phi$')
ylabel('$\dot{m}_\mathrm{NO}/\dot{m}_\mathrm{NO,\mathrm{max}}$')
l=legend(labels,'location','best');
set(l,'interpreter','latex');
hold off

%compute the weight
dt = 200; %[s]
UInf = 200; %[s]
altitude = 43000*ftToM;
ts = 0:dt:100000;
Ws = zeros(1,length(ts));
mDotFs = zeros(1,length(ts));
count=0;
while (HondaJet.airframe.WInitial-HondaJet.airframe.W) < HondaJet.airframe.WFuel
    Ws(count+1)=HondaJet.airframe.W;
    mDotFs(count+1)=fuelMassFlow(UInf, altitude, HondaJet, fluid);
    HondaJet.airframe.W = HondaJet.airframe.W ...
        -mDotFs(count+1)*g*dt;
    count=count+1;
end
ts(count:end)=[];
Ws(count:end)=[];
mDotFs(count:end)=[];
xs=ts*UInf;
range = xs(end);