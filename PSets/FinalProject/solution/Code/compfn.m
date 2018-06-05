function [p_ratio, eta_c, T_out, Work_in] =...
    compfn(cp, gamma, T, eta_st, u_x, U, n, alpha_a, beta_b)

dU_theta = U - u_x * (tan(abs(alpha_a)) + tan(abs(beta_b)));
UdU_theta = U * dU_theta;

T_out = T + n * UdU_theta/cp;
Work_in = n * UdU_theta;

p_ratio = 1;
for i = 1:n
    p_ratio = p_ratio * (1 + eta_st * UdU_theta / ...
        (cp * T + (i - 1) * UdU_theta))^(gamma / (gamma - 1));
end
eta_c = (p_ratio^((gamma - 1)/gamma) - 1) / ...
    (n * UdU_theta/cp/T);

end