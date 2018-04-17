function res = dydt(t, x)
res = zeros(2,1);
res(1) = x(2);
res(2) = exp(-x(1));