% Script: problemSet4Driver
% -------------------------------------------------------------------------
%   This script is used to solve the second problem of problem set 4: the
%   compressor map

clear; clc; close all

%% parameters
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
% pressure ratios
pi.f = 2.0;
% adiabatic efficiencies
n.d = 0.97;
n.f = 0.92;
n.t = 0.91;
n.n = 0.98;
n.n1 = 0.98; %assumed
% angles
alpha.a = 30*3.14159/180;
beta.a = 60*3.14159/180;
alpha.b = 60*3.14159/180;
beta.b = 30*3.14159/180;

%% Freestream conditions
[T0 , a0 , p0 , ~] = atmosisa(h);
U0=a0*M0;

%% ### TODO ###
% Here you need to define and compute:
% p03  : stationary pressure at stage 3 (leaving the compressor)
% n.st : stage efficiency
% rBar : effective mean radius
% T02pt5 : stagnition temperature at stage 2.5
% A02pt5 : area at stage 2.5
% #########################################################################

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