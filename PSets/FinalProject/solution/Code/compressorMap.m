function [pc_ratio, fM2_rev, eta_c, T_out, Work_in] =...
    compressorMap(cp, gamma, T, eta_st, u_x, U, n, alpha_a, beta_b)
R = 287;
[pc_ratio, eta_c, T_out, Work_in] =...
    compfn(cp, gamma, T, eta_st, u_x, U, n, alpha_a, beta_b);

%M025 = u_x/sqrt(gamma * R * T);
M025 = sqrt(1 / (gamma * R * T/u_x^2 - (gamma - 1)/2));
fM2 = (1/M025) * (((gamma + 1)/2)^(-(gamma + 1)/2/(gamma - 1))) * ...
    ((1 + (gamma - 1)/2 * M025^2)^((gamma + 1)/2/(gamma - 1)));
fM2_rev = 1 / fM2;

end