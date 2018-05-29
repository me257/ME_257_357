function [T04, phi, LHV] = Combustor(T03)
% C16H34 + 24.5 (O2 + 3.76 N2) = 16 CO2 + 17 H2O + (O2 + 3.76 N2)
global f gamma
Ru = 8.314e3;
T0 = 298.15;

% molecular weight of each species
W_F = 12 * 16 + 34;
W_O2 = 32;
W_N2 = 28;
W_CO2 = 44;
W_H2O = 18;

% Gas constant [J/kg/K]
R_F = Ru / W_F;
R_O2 = Ru / W_O2;
R_N2 = Ru / W_N2;
R_CO2 = Ru / W_CO2;
R_H2O = Ru / W_H2O;

% equivalence ratio
fst = W_F / (24.5 * (32 + 3.76 * 28));
phi = f / fst;

% stoichiometry coefficients
nu_CO2 = phi * 16;
nu_H2O = phi * 17;
nu_F = phi;
nu_A = 24.5;

% heat of formations [kcal/mol]
hf0_CO2 = -94; 
hf0_H2O = -57.8;
hf0_F = -89.6;

% Cp [J/kg/K]
Cp_F = gamma / (gamma - 1) * R_F;
Cp_O2 = gamma / (gamma - 1) * R_O2;
Cp_N2 = gamma / (gamma - 1) * R_N2;
Cp_CO2 = gamma / (gamma - 1) * R_CO2;
Cp_H2O = gamma / (gamma - 1) * R_H2O;
W_tot = nu_F * W_F + nu_A * (W_O2 + 3.76 * W_N2);
Yr_F = nu_F * W_F / W_tot;
Yr_O2 = nu_A * W_O2 / W_tot;
Yr_N2 = 3.76 * nu_A * W_N2 / W_tot;
Yp_CO2 = nu_CO2 * W_CO2 / W_tot;
Yp_H2O = nu_H2O * W_H2O / W_tot;
Yp_O2 = nu_A * (1 - phi) * W_O2 / W_tot;
Yp_N2 = 3.76 * nu_A * (1 - phi) * W_N2 / W_tot;
Cp_r = Yr_F * Cp_F + Yr_O2 * Cp_O2 + Yr_N2 * Cp_N2;
Cp_p = Yp_CO2 * Cp_CO2 + Yp_H2O + Cp_H2O + Yp_O2 * Cp_O2 + Yp_N2 * Cp_N2;

kcal2J = 4184;
LHV = -(nu_CO2 * hf0_CO2 + nu_H2O * hf0_H2O - nu_F * hf0_F);
Hf = kcal2J * LHV / (W_F * 1e-3);

T04 = T0 + (Cp_r * (T03 - T0) + Hf) / Cp_p;