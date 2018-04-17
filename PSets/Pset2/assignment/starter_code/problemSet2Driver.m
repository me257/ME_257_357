clear; clc; close all
%% problem 2
%define a structure containing the HondaJet's parameters
HondaJet.engine.overallCompressionRatio = 24;   %[]
HondaJet.engine.turbineInletTemperature=1600;   %[K]
HondaJet.engine.mDotA=1.7; 
HondaJet.engine.combustor.betaMix=100e-6;       %[]
HondaJet.engine.combustor.mDotA=1.27;           %[kg/s]
HondaJet.engine.combustor.CS=0.2;               %[]
HondaJet.engine.combustor.rF=810;               %[kg/m^3]
HondaJet.engine.combustor.D0=50e-6;             %[]
HondaJet.engine.combustor.beta=0.5e-6;          %[m^2]
HondaJet.engine.combustor.nDil=25;              %[]
HondaJet.engine.combustor.dDil=1e-2;            %[m]
HondaJet.engine.combustor.CDil=5e-3;            %[]
HondaJet.engine.combustor.Ta=2100;              %[K]
HondaJet.engine.combustor.n=1.0;                %[]
HondaJet.engine.combustor.m=0.2;                %[]
HondaJet.engine.combustor.A=240;                %[Pa.s]
HondaJet.engine.combustor.Ro=0.25;              %[m]
HondaJet.engine.combustor.Ri=0.1;               %[m]
HondaJet.engine.combustor.Ue=0.75;            %[m/s]

%Set up the fluid
fluid.gas = Solution('ndodecane_mech.cti');
fluid.cp = 1005; %[J/kg]
fluid.R = 287; %[J/kg.K]
iC12H26 = 1;
iO2 = 2;
iN2 = 3;
fluid.LHV = Hc(phi==1)*1000/W(iC12H26); %[J/kg]
mC12H26st = mass.reactants(iC12H26,phi==1);
mO2st = mass.reactants(iO2,phi==1);
mN2st = mass.reactants(iN2,phi==1);
fluid.fst = mC12H26st/(mO2st+mN2st);

%compute the length scales
HondaJet.engine.combustor.betaC =linspace(1,10);
L = combustorDesign(200, 43000*ftToM, HondaJet.engine, fluid);

%Make the plots
figure()
semilogy(HondaJet.engine.combustor.betaC,L.LTotal)
hold on
semilogy(HondaJet.engine.combustor.betaC,L.LEvap)
semilogy(HondaJet.engine.combustor.betaC,L.LDil)
semilogy(HondaJet.engine.combustor.betaC,L.LSrz)
hold off
set(gca,'FontSize',16)
labels = {'$L_\mathrm{total}$'
          '$L_\mathrm{evap}$'
          '$L_\mathrm{dil}$'
          '$L_\mathrm{srz}$'};
xlabel('$\beta_\mathrm{c}$','FontSize',16)
ylabel('$L\ [\mathrm{m}]$','FontSize',16)
leg = legend(labels,'location','best');
set(leg,'interpreter','latex','FontSize',16)
xlim([1 10])
