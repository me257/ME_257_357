%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ME357 Spring 2017
%  Final Project
%  3. Turbine Map 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Given properties
cp = 1005;
gamma = 1.4;
T = 1600;
ux_ref = 732;
U_ref = 216;
n = 1;
alpha_b = 60 / 180 * pi;
beta_c  = 45 / 180 * pi;
eta_st = @(ux,uu) (0.95 - 0.03 * (ux/uu - ux_ref/U_ref).^2 -...
    0.2 * abs(uu - U_ref)/U_ref);
r_mean = 0.08;

% Preparation for the turbine map
N = 15000:1000:30000;
ux = 200:10:750;
lN = length(N);
lux = length(ux);
map_mDot = zeros(lux,1);
map_pRatio = zeros(lux,1);
map_etat = zeros(lux,1);

% Calculate the turbine map
for j = 1:lN
    U = N(j)/60 * 2 * pi * r_mean;
    for i = 1:lux
        [p_ratio, fM4_rev, eta_t, T_out, Work_out] =...
            turbineMap(cp, gamma, T, eta_st(ux(i),U), ux(i), U, n, alpha_b, beta_c);
        map_mDot(i) = fM4_rev;
        map_pRatio(i) = 1 / p_ratio;
        map_etat(i) = eta_t;

    end

    figure(1)
    plot(map_mDot,map_pRatio, 'k');
    hold on;
    
    figure(2)
    plot(map_mDot * N(j),map_pRatio, 'k');
    hold on;
    mDot_contour(j,:) = map_mDot * N(j);
    pRatio_contour(j,:) = map_pRatio;
    etat_contour(j,:) = map_etat;
end

% Plot turbine map
figure(1)
xlabel('$1/f(M_{04})$','FontSize',22,'FontWeight','bold','Interpreter','latex');
ylabel('$p_{04}/p_{04.5}$','FontSize',22,'FontWeight','bold','Interpreter','latex');
axesh = findobj('Type', 'axes');
set(axesh, 'Box','on');
set(gca,'FontSize',20);

figure(2)
eta_values = 0.5:0.01:0.98;
contour(mDot_contour, pRatio_contour, etat_contour, eta_values, 'LineWidth',3.0);
cbar = colorbar;
xlabel('$N * 1/f(M_{04})$','FontSize',22,'FontWeight','bold','Interpreter','latex');
ylabel('$p_{04}/p_{04.5}$','FontSize',22,'FontWeight','bold','Interpreter','latex');
ylabel(cbar,'Turbine Efficiency, $\eta_t$',...
    'FontSize',22,'FontWeight','bold','Interpreter','latex');
axesh = findobj('Type', 'axes');
set(axesh, 'Box','on');
set(gca,'FontSize',20);