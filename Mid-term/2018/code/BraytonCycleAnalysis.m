clear; close all; clc
%% Constants and known parameters
global Cp gamma
gamma = 1.4;
R = 285;
Cp = gamma * R / (gamma - 1);

M0 = 0.85;
T0 = 220;
p0 = 24e3;

% core air flow rate
mdot_ac = 3.0;
% bypass ratio
beta = 2.9;
% fuel-to-air ratio
global f
f = 0.0043;
% number of compressor stages
nstage = 3;
% rotational velocity
global Omega
Omega = 20000 * 6 / 5; % rpm
% vortex design constant
global b
b = 175 * 5 / 6;

p025p02 = 2.0;

eta_d = 0.95;
eta_f = 0.92;
eta_c = 0.87;
eta_t = 0.91;
eta_n = 0.98;

%% Pressure and temperature ratios
T0bT0a_inc = @(x, eta) 1 + (1/eta) * (x^((gamma - 1)/gamma) - 1);
p0bp0a_inc = @(x, eta) (1 + eta * (x - 1))^(gamma / (gamma - 1));
T0bT0a_dec = @(x, eta) 1 - eta * (1 - x^((gamma - 1)/gamma));
p0bp0a_dec = @(x, eta) (1 - (1/eta) * (1 - x))^(gamma/(gamma - 1));

%% Brayton cycle analysis
% Diffuser
T02T0 = 1 + (gamma - 1) / 2 * M0^2;
p02p0 = p0bp0a_inc(T02T0, eta_d);
T02 = T02T0 * T0;
p02 = p02p0 * p0;

% Fan
T025T02 = T0bT0a_inc(p025p02, eta_f);
T025 = T025T02 * T02;
p025 = p025p02 * p02;

% High-Pressure Compressor
p03p025_ = HighPressorCompressor(T025);
p03p025 = p03p025_^nstage; % 24; % 
T03T025 = T0bT0a_inc(p03p025, eta_c);
T03 = T03T025 * T025;
p03 = p03p025 * p025;

% Combustor
[T04, phi, LHV] = Combustor(T03);
% T04 = 1600;
p04 = p03;
T04T03 = T04 / T03;
p04p03 = p04 / p03;

% High-Pressure Turbine
T045T04 = 1 - 1 / (1 + f) * (T03 - T025) / T04;
p045p04 = p0bp0a_dec(T045T04, eta_t);
T045 = T045T04 * T04;
p045 = p045p04 * p04;

% Low-Pressure Turbine
T05T045 = 1 - beta / (1 + f) * (T025 - T02) / T045;
p05p045 = p0bp0a_dec(T05T045, eta_t);
T05 = T05T045 * T045;
p05 = p05p045 * p045;

% Core Nozzle
p08p05 = p0 / p05;
T08T05 = T0bT0a_dec(p08p05, eta_n);
p08 = p08p05 * p05;
T08 = T08T05 * T05;

% Fan/Bypass Nozzle
p018p02 = p0 / p02;
T018T02 = T0bT0a_dec(p018p02, eta_n);
p018 = p018p02 * p02;
T018 = T018T02 * T02;

%% Engine parameters
U1e = sqrt(2 * eta_f * Cp * T025 * (1 - (p0 / p025)^((gamma - 1)/gamma)));
Ue  = sqrt(2 * eta_n * Cp * T05  * (1 - (p0 / p05 )^((gamma - 1)/gamma)));

U0 = M0 * sqrt(gamma * R * T0);
T = mdot_ac * ((1 + f) * Ue + beta * U1e - (1 + beta) * U0);

TSFC = f / ((1 + f) * Ue + beta * U1e - (1 + beta) * U0);

eta_p = (2 * U0 * ((1 + f) * Ue + beta * U1e - (1 + beta) * U0))/ ...
    ((1 + f) * Ue^2 + beta * U1e^2 - (1 + beta) * U0^2);