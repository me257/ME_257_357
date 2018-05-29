function p03p02 = HighPressorCompressor(T01)
global gamma Cp Omega b
eta_st = 0.88;
p03p02 = (1 + eta_st * (Omega/60 * 2 * b) / (Cp * T01))^(gamma / (gamma - 1));