function [p_ratio, eta_t, T_out, Work_out] = ...
    turbfn(cp, gamma, T, eta_st, u_x, U, n, alpha_b, beta_c)

dU_theta = u_x * (tan(abs(alpha_b)) + tan(abs(beta_c))) - U;
UdU_theta = U * dU_theta;

T_out = T - n * UdU_theta/cp;
Work_out = n * UdU_theta;

p_ratio = 1;
for i = 1:1:n
    p_ratio = p_ratio * ((1 - 1/eta_st * UdU_theta/...
        (cp * T - (i - 1) * UdU_theta))^(gamma/(gamma-1)));
end

eta_t = (n * UdU_theta/cp/T) / (1 - p_ratio^((gamma - 1)/gamma));

end