clear; clc; close all

%% problem 1.1b
ftToM = 12*0.0254;
phi = 0:0.1:2;
%create a matrix of the moles of each species (nSp X nPhi)
nu.reactants =         [phi                            %C12H26 [mol]
                        37/2*ones(1,length(phi))        %O2 [mol]
                        37/2*3.76*ones(1,length(phi))   %N2 [mol]
                        zeros(1,length(phi))            %CO2[mol]
                        zeros(1,length(phi))];          %H2O[mol]
nu.products =          [max(0,phi-1)                    %C12H26 [mol]
                        37/2*max(0,1-phi)               %O2 [mol]
                        37/2*3.76*ones(1,length(phi))   %N2 [mol]
                        12*min(1,phi)                   %CO2[mol]
                        13*min(1,phi)];                 %H2O[mol]        
%molecular weights
W = [12*12+26*1     %C12H26 [g/mol]
    2*16            %O2 [g/mol]
    2*14            %N2 [g/mol]
    12+16*2         %CO2 [g/mol]
    1*2+16];        %H2O [g/mol]
%compute the masses
iCO2 = 4; iH2O=5;
mass.reactants = (W*ones(1,length(phi))).*nu.reactants; %[g]
mass.products =  (W*ones(1,length(phi))).*nu.products; %[g]

%make the plots
figure()
plot(phi,mass.products(iCO2,:)./sum(mass.products,1))
hold on
plot(phi,mass.products(iH2O,:)./sum(mass.products,1))
hold off
xlabel('$\phi$')
ylabel('$Y$')
l=legend('$Y_\mathrm{CO_2}$','$Y_\mathrm{H_2O}$','location','best');
set(l,'interpreter','latex')

%% problem 1.2a
%heat of combustion
Hf = [-1.71e3       %C12H26 [J/g]
     0              %O2 [J/g]
     0              %N2 [J/g]
     -8.94e3        %CO2 [J/g]
     -13.4e3];      %H2O [J/g]
Hc = Hf'*(mass.reactants-mass.products); %[J]
%heat capacities
cp = [3.37          %C12H26 [J/g.K]
     1.09           %O2 [J/g.K]
     1.18           %N2 [J/g.K]
     1.22           %CO2 [J/g.K]
     2.38];         %H2O [J/g.K]
Cp.reactants = cp'*mass.reactants; %[J/K]
Cp.products = cp'*mass.products; %[J/K]
%flame temperature
TReactants = 930; %[K]
TRef = 298; %[K]
TAdiabatic = TRef+(Hc+Cp.reactants*(TReactants-TRef))./Cp.products; %[K]
%compare using cantera
pReactants = 5.7e5; %[bar]
nPhi = length(phi);
gas = Solution('ndodecane_mech.cti');
nSp = nSpecies(gas);
canteraIndices =   [speciesIndex(gas,'nc12h26')
                    speciesIndex(gas,'o2')
                    speciesIndex(gas,'n2')
                    speciesIndex(gas,'co2')
                    speciesIndex(gas,'h2o')];
TCantera = zeros(1,nPhi);
for iPhi=1:nPhi
    YReactants = zeros(nSp,1);
    YReactants(canteraIndices) = mass.reactants(:,iPhi);
    set(gas,'T',TReactants,'P',pReactants,'Y',YReactants);
    equilibrate(gas,'HP');
    TCantera(iPhi)=temperature(gas);
end
% plotting
figure()
plot(phi,TAdiabatic)
hold on
plot(phi,TCantera)
hold off
xlabel('$\phi$')
ylabel('$T_\mathrm{ad}\ \mathrm{[K]}$')
l=legend('calorically perfect','cantera','location','best');
set(l,'interpreter','latex')

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
