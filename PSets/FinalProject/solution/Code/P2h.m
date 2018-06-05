%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ME357 Spring 2017
%  Final Project
%  2. Turbine Design (h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Given properties
cp = 1005;
gamma = 1.4;
T = 1600;
eta_st = 0.95;
u_x = 731.94;
U = 215.65;
n = 1;

% Turbine geometries
alpha_bd = 20:1:85;
beta_cd = 20:1:85;
alpha_b = alpha_bd / 180 * pi;
beta_c  = beta_cd / 180 * pi;
lb = length(alpha_b);
lc = length(beta_c); 

% Calculate contours
work = zeros(lc,lb);
for i = 1:lc
    for j = 1:lb
        [p_ratio, eta_t, T_out, Work_out] =...
            turbfn(cp, gamma, T, eta_st, u_x, U, n, alpha_b(j), beta_c(i));
        work(i,j) = Work_out;
    end
end
work_values = logspace(4,7,20);
work_values = [2.71e5 work_values];

% Plot contours
figure(1)
co = [1 0 0;0 1 0;0 0 1];
contour(alpha_bd, beta_cd, work, work_values,'ShowText','on','LineWidth',2.0)
hold on
colorbar
xlabel('$\alpha_b [^o]$','FontSize',25,'FontWeight','bold','Interpreter','latex');
ylabel('$\beta_c [^o]$','FontSize',25,'FontWeight','bold','Interpreter','latex');
axesh = findobj('Type', 'axes');
set(axesh, 'Box','on');
set(gca,'FontSize',23);
plot(60, 45, 'ro','MarkerSize',10.0, 'LineWidth',2.0)
text(61,45,'Turbine Design Point','Color','r','FontSize',15.0, ...
    'FontWeight','Bold','Interpreter','latex')
