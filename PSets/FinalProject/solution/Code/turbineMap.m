function [pt_ratio, fM4_rev, eta_t, T_out, Work_out] =...
    turbineMap(cp, gamma, T, eta_st, u_x, U, n, alpha_b, beta_c)

R = 287;
[pt_ratio, eta_t, T_out, Work_out] =...
    turbfn(cp, gamma, T, eta_st, u_x, U, n, alpha_b, beta_c);

M4 = sqrt(1 / (gamma * R * T/u_x^2 - (gamma - 1)/2));
fM4 = (1/M4) * (((gamma + 1)/2)^(-(gamma + 1)/2/(gamma - 1))) * ...
    ((1 + (gamma - 1)/2 * M4^2)^((gamma + 1)/2/(gamma - 1)));
fM4_rev = 1 / fM4;

end