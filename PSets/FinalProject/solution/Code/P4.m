%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ME357 Spring 2017
%  Final Project
%  4. Compressor-Turbine Matching and Operating Line
%
%  Rui Xu (ruixu@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Fluid properties
cp = 1005;
gamma = 1.4;
R = 287;
LHV = 43e6;

%% Compressor Map
% Given for compressor
n_c = 4;
T025 = 311;
p025 = 82.4e3;
ux_c_ref = 75;
U_c_ref = 310;
mDot_c = 2.825;
r_mean_c = 0.115;
eta_c_st = @(x) (0.87 - 16 * (x - ux_c_ref/U_c_ref).^2);
alpha_a_comp = 30 / 180 * pi;
beta_b_comp  = 30 / 180 * pi;

% Preparation for compressor map
N = 15000:1000:30000;
ux_c = 20:1:130;
lN = length(N);
lux_c = length(ux_c);
Cmap_mDot = zeros(lux_c,1);
Cmap_pRatio = zeros(lux_c,1);
Cmap_etac = zeros(lux_c,1);

% Calculate and plot the compressor map
figure(1)
hold on
for j = 1:lN
    U_c = N(j)/60 * 2 * pi * r_mean_c;
    for i = 1:lux_c
        [pc_ratio, fM2_rev, eta_c, T_out, Work_in] =...
            compressorMap(cp, gamma, T025, eta_c_st(ux_c(i)/U_c),...
            ux_c(i), U_c, n_c, alpha_a_comp, beta_b_comp);
        Cmap_mDot(i) = fM2_rev;
        Cmap_pRatio(i) = pc_ratio;
        Cmap_etac(i) = eta_c;
    end
    
    % Cut off at surge line
    Cmap_etac = Cmap_etac((find(Cmap_pRatio == max(Cmap_pRatio))):end);
    Cmap_mDot = Cmap_mDot((find(Cmap_pRatio == max(Cmap_pRatio))):end);
    Cmap_pRatio = Cmap_pRatio((find(Cmap_pRatio == max(Cmap_pRatio))):end);
    plot(Cmap_mDot,Cmap_pRatio, 'k');
    
    % Prepare for compressor efficiency contour plot
    mDot_c_contour(j,:) = linspace(Cmap_mDot(1), Cmap_mDot(end), 100);
    pRatio_c_contour(j,:) = interp1( Cmap_mDot,Cmap_pRatio, mDot_c_contour(j,:) );
    etac_contour(j,:) = interp1( Cmap_mDot,Cmap_etac, mDot_c_contour(j,:) );
end

% Plot efficiency contour
figure(1)
eta_values = 0.1:0.01:0.9;
contour(mDot_c_contour, pRatio_c_contour, etac_contour, eta_values, 'LineWidth',3.0)
cbar = colorbar;
scale = axis;
axis([0.1 0.7 0 12]);
xlabel('$1/f(M_{02.5})$','FontSize',22,'FontWeight','bold','Interpreter','latex');
ylabel('$p_{03}/p_{02.5}$','FontSize',22,'FontWeight','bold','Interpreter','latex');
ylabel(cbar,'Compressor Efficiency, $\eta_c$',...
    'FontSize',22,'FontWeight','bold','Interpreter','latex');
axesh = findobj('Type', 'axes');
set(axesh, 'Box','on');
set(gca,'FontSize',20);

%% Turbine Map
% Given for turbine
T04_t = 1600;
ux_t_ref = 732;
U_t_ref = 216;
n_t = 1;
alpha_b_turb = 60 / 180 * pi;
beta_c_turb  = 45 / 180 * pi;
eta_t_st = @(ux,uu) (0.95 - 0.03 * (ux/uu - ux_t_ref/U_t_ref).^2 -...
    0.2 * abs(uu - U_t_ref)/U_t_ref);
r_mean_t = 0.08;

% Preparation for turbine map
ux_t = 200:10:750;
lux_t = length(ux_t);
Tmap_mDot = zeros(lux_t,1);
Tmap_pRatio = zeros(lux_t,1);
Tmap_etat = zeros(lux_t,1);

% Calculate and plot turbine map
for j = 1:lN
    U_t = N(j)/60 * 2 * pi * r_mean_t;
    for i = 1:lux_t
        [pt_ratio, fM4_rev, eta_t, T_out, Work_out] =...
            turbineMap(cp, gamma, T04_t, eta_t_st(ux_t(i),U_t), ...
            ux_t(i), U_t, n_t, alpha_b_turb, beta_c_turb);
        Tmap_mDot(i) = fM4_rev;
        Tmap_pRatio(i) = 1 / pt_ratio;
        Tmap_etat(i) = eta_t;
    end
    
    figure(2)
    hold on;
    plot(Tmap_mDot * N(j),Tmap_pRatio, 'k');
    mDot_t_contour(j,:) = Tmap_mDot * N(j);
    pRatio_t_contour(j,:) = Tmap_pRatio;
    etat_contour(j,:) = Tmap_etat;
end

% Plot efficiency contour
figure(2)
eta_values = 0.5:0.01:0.98;
contour(mDot_t_contour, pRatio_t_contour, etat_contour, eta_values, 'LineWidth',3.0);
cbar = colorbar;
xlabel('$N * 1/f(M_{04})$','FontSize',22,'FontWeight','bold','Interpreter','latex');
ylabel('$p_{04}/p_{04.5}$','FontSize',22,'FontWeight','bold','Interpreter','latex');
ylabel(cbar,'Turbine Efficiency, $\eta_t$',...
    'FontSize',22,'FontWeight','bold','Interpreter','latex');
axesh = findobj('Type', 'axes');
set(axesh, 'Box','on');
set(gca,'FontSize',20);

%% Operating line
% Areas directly calculated from compressor and turbine geometries.
A025 = 0.04154756284;
A04 = 0.002826930733;

% Preparation for operating line calculation
Op_c_mDot = zeros(1, lN);
Op_c_pRatio = zeros(1, lN);
Op_t_mDot = zeros(1, lN);
Op_t_pRatio = zeros(1, lN);
for j = 1:lN
    U_c = N(j)/60 * 2 * pi * r_mean_c;
    U_t = N(j)/60 * 2 * pi * r_mean_t;
    fun = 1e6;
    tol = 0.01;
    
    % Guess a ux in compressor
    for ux_c_guess = 20:0.1:300
        
        % Get compressor map output
        [pc_ratio, fM025_rev, eta_c, T03_guess, Work_in] =...
            compressorMap(cp, gamma, T025, eta_c_st(ux_c_guess/U_c),...
            ux_c_guess, U_c, n_c, alpha_a_comp, beta_b_comp);
        p03_guess = pc_ratio * p025;
        p04_guess = p03_guess;
        mDot_c_guess = gamma * (((gamma + 1)/2)^(-(gamma + 1)/2/(gamma - 1))) * ...
            p025 * A025 * fM025_rev/sqrt(gamma * R * T025);
        
        % work balance
        wc = Work_in;
        wt = wc;
        ux_t_guess = (wt/U_t + U_t) / ...
            (tan(abs(alpha_b_turb)) + tan(abs(beta_c_turb)));
        T04_guess = ux_t_guess^2/R * (gamma + 1)/2/gamma;
        
        % Get turbine map output
        [pt_ratio, fM04_rev, eta_t, T05_guess, Work_out] =...
            turbineMap(cp, gamma, T04_guess, eta_t_st(ux_t_guess,U_t), ...
            ux_t_guess, U_t, n_t, alpha_b_turb, beta_c_turb);
        f_guess = (T04_guess/T03_guess - 1) / ...
            (LHV/(cp * T03_guess) - T04_guess/T03_guess);
        mDot_t_guess = gamma * (((gamma + 1)/2)^(-(gamma + 1)/2/(gamma - 1))) * ...
            p04_guess * A04 * fM04_rev/sqrt(gamma * R * T04_guess);
        
        % Check mass balance
        fun = mDot_t_guess - mDot_c_guess * (1 + f_guess);
        if (abs(fun) < tol)
            break;
        end
    end
    
    Op_c_mDot(j) = fM025_rev;
    Op_c_pRatio(j) = pc_ratio;
    Op_t_mDot(j) = fM04_rev;
    Op_t_pRatio(j) = 1 / pt_ratio;
end

cc = [0.3091,1.66;
    0.325,1.67;
    0.3505,1.68;
    0.37,1.7;
    0.4,1.8;
    0.5,2.1;
    0.55,2.5;
    0.5934,3.5];

tt = [0.9011,1.25;
    1,1.26;
    1.25,1.335;
    1.5,1.423;
    1.75,1.547;
    2,1.687;
    2.25,1.85;
    2.5,2.1;
    2.7,2.3];

xxc = linspace(0.3091,0.5934,16);
xxt = linspace(0.9011,2.7,16);
yyc = interp1(cc(:,1),cc(:,2),xxc,'spline');
yyt = interp1(tt(:,1),tt(:,2),xxt,'spline');

figure(1)
plot(xxc,yyc,'r--','LineWidth',3.0)

figure(2)
plot(xxt * 1e4, yyt,'r--','LineWidth',3.0)