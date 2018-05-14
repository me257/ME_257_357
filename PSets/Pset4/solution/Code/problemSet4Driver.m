% Script: problemSet4Driver
% -------------------------------------------------------------------------
%   This script is used to solve the second problem of problem set 4: the
%   compressor map

clear; clc; close all

%parameters
Thrust=3400; %[N]
BPR = 2.9;
T04 = 1600; %[k]
M0 = 0.7;
gamma = 1.4; %assumed
cp = 1005; %assumed [J/kg.K]
R = 8314/28.85; %assumed [J/kg.K]
uxRef = 75; %[m/s]
URef = 307; %[m/s]
h = 30000*12*0.0254; %[m]
nStages = 4;
rHubTipRatio = 0.6;
%pressure ratios
pi.f = 2.0;
%adiabatic efficiencies
n.d = 0.97;
n.f = 0.92;
n.t = 0.91;
n.n = 0.98;
n.n1 = 0.98; %assumed
%angles
alpha.a = 30*3.14159/180;
beta.a = 60*3.14159/180;
alpha.b = 60*3.14159/180;
beta.b = 30*3.14159/180;

%Freestream conditions
[T0 , a0 , p0 , ~] = atmosisa(h);
U0=a0*M0;

%get turbine inlet temperature ratio
tau.l = T04./T0;

%compute the reference conditions
tau.r = 1+(gamma-1)./2.*M0.^2;
pi.r = tau.r.^(gamma/(gamma-1.0));

%diffuser pressure ratio
tau.d = 1.0; %by assumption
pi.d = (1.0+n.d.*(tau.r-1.0)).^(gamma./(gamma-1.0))./pi.r;

%compute the fan temperature ratio
tau.f = 1.0+1.0./n.f.*(pi.f.^((gamma-1.0)./gamma)-1.0);
%compute compressor efficiency and compressor temperature ratio
uxoU = uxRef/URef;
n.st = @(uxoU) 0.87-16*(uxoU-uxRef/URef).^2;
%determine the temperature ratio
T02pt5=T0*tau.r*tau.d*tau.f;
p02pt5 = p0*pi.r*pi.d*pi.f;
T03=T02pt5;
p03=p02pt5;
for nStage = 1:nStages
    Duth = URef*(1.0-uxoU*(tan(alpha.a)+tan(beta.b)));
    T03 = T03*(1.0+URef*Duth/(cp*T03));
    p03 = p03*(1.0+n.st(uxoU)*URef*Duth/(cp*T03))^(gamma/(gamma-1));
end
tau.c=T03/T02pt5;
pi.c=p03/p02pt5;

%compute the turbine temperature ratio and pressure ratio
tau.hpt = 1.0+tau.r.*tau.d.*tau.f.*(1.0-tau.c)./tau.l;
tau.lpt = 1.0+tau.r.*tau.d.*(1-tau.f).*(1.0+BPR)./(tau.l.*tau.hpt); %lowpressure turbine
tau.t = tau.hpt.*tau.lpt;
pi.t = (1.0+(n.t.*(tau.t-1.0))).^(gamma./(gamma-1));

%calculate exit velocities
Ue = a0.*sqrt(2.*n.n.*tau.l.*tau.t./(gamma-1.0).*(1.0-(pi.r.*pi.d.*pi.f.*pi.c.*pi.t).^((1-gamma)./gamma)));
Ue1 = a0.*sqrt(2.0.*n.n1.*tau.r.*tau.f./(gamma-1.0).*(1.0-(pi.r.*pi.d.*pi.f).^((1.0-gamma)./gamma)));

%calculate the mass flow
mDota =  Thrust/(Ue-U0+BPR.*(Ue1-U0));

%caluculate the area
GStar = @(T0,p0) 0.04042*p0/sqrt(T0); %mass flux at the choked condition for gamma=1.4
AoAstr = @(M) 1/M*(2/(gamma+1)*(1+(gamma-1)/2*M^2)).^((gamma+1)/2/(gamma-1)); %ratio of area to that at the choked condtino
T2pt5 = T02pt5-uxRef^2/(2*cp); %subtract the kinetic energy from the stagnation temperature
M2pt5 = uxRef/(sqrt(gamma*R*T2pt5));
A2pt5 = mDota/GStar(T02pt5,p02pt5)*AoAstr(M2pt5); %area at station 2pt5
%calculate the mean radius
rTip = sqrt(A2pt5/(1-rHubTipRatio^2)/3.14159);
rHub = rTip*rHubTipRatio;
rBar = (rTip+rHub)/2;

%determine the stagnation pressure out of the compressor
p03 = p02pt5*pi.c;

%% Compressor Map

%compute the parameters needed for the compressor map function
args.rBar=rBar;
args.URef=URef;
args.gamma=gamma;
args.R = R;
args.cp = cp;
args.T02pt5 = T02pt5;
args.A2pt5 = A2pt5;
args.n = n;
args.alpha = alpha;
args.beta = beta;
args.nStages = nStages;
args.uxRef = uxRef;
args.T04 = T04;

%call an abstracted function to do the plotting for the compressor map
plotCompressorMap(args)