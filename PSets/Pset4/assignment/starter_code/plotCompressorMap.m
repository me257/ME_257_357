function plotCompressorMap(args)
    % function: plotCompressorMap
    % -------------------------------------------------------------------------
    %   This function is used to plot the compressor map    

    %unpack the arguments required for this function
    rBar=args.rBar;
    URef=args.URef;
    gamma=args.gamma;
    R=args.R;
    cp = args.cp;
    T02pt5 = args.T02pt5;
    A2pt5 = args.A2pt5;
    n = args.n;
    alpha = args.alpha;
    beta = args.beta;
    nStages = args.nStages;
    uxRef = args.uxRef;
    T04 = args.T04;
    
    %compute the rotation velocity
    N = URef/(2*3.14159*rBar)*60; %[RPM]

    %define the stage 
    Ntemp = [N , (15000: 1000: 30000) ];
    Ntemp = sort(Ntemp);
    uxrange = linspace(20,130,101);
    x = zeros(size(uxrange));
    y = zeros(size(uxrange));
    eta = zeros(size(uxrange));
    xmax = zeros(size(Ntemp));
    ymax = xmax;
    uxsurge = xmax;

    eff_nx = 75;
    eff_ny = 50;

    N_plot_interval = max(1, floor(length(Ntemp)/6));
    surge_plot_coeffs = [;];

    figure(1);
    for j = 1:length(Ntemp)
        Utemp = Ntemp(j)/60*2*pi*rBar;
        for i = 1:length(uxrange)
            temp = uxrange(i);
            [y(i), x(i), eta(i)] = ...
                compmap(gamma,R,cp,T02pt5,A2pt5,n.st,temp,Utemp,nStages,alpha.a,beta.b);
        end

        %find the maximum point
        eta = eta((find(y==max(y))):end);
        x = x((find(y==max(y))):end);
        y = y((find(y==max(y))):end);

        if length(x) > 1
            x_matrix(j,:) = linspace(x(1), x(end), eff_nx);
            y_matrix(j,:) = interp1( x,y, x_matrix(j,:) );
            eta_matrix(j,:) = interp1( x,eta, x_matrix(j,:) );
        else
            x_matrix(j,:) = repmat(x, 1, eff_nx);
            y_matrix(j,:) = repmat(y, 1, eff_nx);
            eta_matrix(j,:) = repmat(eta, 1, eff_nx);
        end


        xmax(j) = x(1);
        ymax(j) = y(1);
        xmax(j) = x((find(y==max(y))));
        ymax(j) = y((find(y==max(y))));

        if ( mod(j-1,N_plot_interval) == 0 || Ntemp(j) == N || j == length(Ntemp) )
        plot(x,y, 'k')
        text(x(end),y(end),['$\leftarrow N/(T_{02.5})^{1/2} = $', num2str(Ntemp(j)/sqrt(T02pt5))],...
         'HorizontalAlignment','left', 'FontWeight', 'Bold', 'FontSize', 14)
     
        if length(x) > 1
            surge_plot_coeffs = [surge_plot_coeffs, [xmax(j); ymax(j)] ];
        end

        end

        hold on
    end

    etalevels = [ (0:0.1:0.5) , (0.55:0.05:0.7) , (0.75:0.01:0.9) ];
    contour(x_matrix, y_matrix, eta_matrix, etalevels, 'LineWidth', 2)
    colormap(winter)
    colorbar


    %curvefit the surge pressure for later use in the rootfinder
    uxcoff = polyfit(Ntemp,ymax,1);
    axis([0.0 0.9 0 1.1*max(ymax)])
    plot(surge_plot_coeffs(1,:) , surge_plot_coeffs(2,:) , ':k', 'LineWidth', 4)
    text(xmax(end-1),ymax(end-1),'Surge line ',...
         'HorizontalAlignment','right', 'FontWeight', 'Bold', 'FontSize', 14)
    [yA, xA] =...
        compmap(gamma,R,cp,T02pt5,A2pt5,n.st,uxRef,URef,nStages,alpha.a,beta.b);
    plot(xA,yA,'o', 'MarkerSize', 10)
    text(xA,yA,'   Point A',...
         'HorizontalAlignment','left', 'FontWeight', 'Bold', 'FontSize', 14)
     a = yA/xA;
     arange = [0.4*a 0.6*a 0.8*a a 1.1*a];
     x = linspace(0.2,0.6,10);
     for i=1:length(arange)
         plot(x,arange(i)*x, 'Color', [0.75 0.075 0.075]);
         text(x(1),x(1)*arange(i),[  '$T_{04}/T_{02.5} =$' num2str((arange(i)/a)*T04/T02pt5)],...
         'HorizontalAlignment','right', 'FontWeight', 'Bold', 'FontSize', 14)
     end
     ylabel('$P_{03}/P_{02.5}$')
     xlabel('$\dot{m}_c  (T_{02.5})^{1/2}/(0.04 p_{02.5} A_{2.5})$')
     title('HPC Compressor Map - Mach 0.7, 30,000 ft')


    %Finding thrust at some operating point

    %Cruise
    ux_cruise = 75;
    T4_target = 1000;
    Thrust_Required = 1200;
    [yB, xB, etac, N_target] = operating_point_N_match_drag(9144, 0.7, Thrust_Required , T4_target, gamma,R,cp,A2pt5,n.st,ux_cruise,nStages,alpha.a,beta.b,rBar)

    U_cruise = N_target /60*2*pi*rBar;
    [yBB xBB] = operating_point(9144, 0.7, Thrust_Required , T4_target, gamma,R,cp,A2pt5,n.st,ux_cruise,U_cruise,nStages,alpha.a,beta.b)

    disp('Found Cruise condition...')

    plot(xB,yB,'rd', 'MarkerSize', 10)
    text(xB,yB,'CRUISE   ',...
         'HorizontalAlignment','right', 'FontWeight', 'Bold', 'FontSize', 14)


    %Now do Sea-level Static Thrust
    ux_slst = 40;
    U_slst = 20000 /60*2*pi*rBar;    % URef;
    [yC, xC] = operating_point(0, 0, 900 , 1600, gamma,R,cp,A2pt5,n.st,ux_slst,U_slst,nStages,alpha.a,beta.b);

    plot(xC,yC,'ms', 'MarkerSize', 10)
    text(xC,yC,'SLST   ',...
         'HorizontalAlignment','right', 'FontWeight', 'Bold', 'FontSize', 14)

end

function [y,x, ec] =compmap(gam,R,cp,T,A2,etast,uxRef,URef,nStages,alfa1,beta2)
    M2 = uxRef/sqrt(gam*R*T);
    Aratio = 1/M2*(2/(gam+1)*(1+(gam-1)/2*M2^2)).^((gam+1)/2/(gam-1));
    Astar = A2/Aratio;
    x = 1/Aratio;
    eta = feval(etast,uxRef/URef);
    [y, ec] = compfn(cp,T,eta,uxRef,URef,nStages,alfa1,beta2);
end

function [pratio, eta_compressor] = compfn(cp,T,eta,uxRef,URef,nStages,alfa1,beta2)

U_c_theta = URef^2/cp/T*(1-uxRef/URef*(tan(alfa1)+tan(abs(beta2))));
pratio = ( 1 + (eta*U_c_theta) ) ^ (3.5*nStages);
eta_compressor = ( (1 + (eta*U_c_theta))^nStages - 1) / ( (1 + U_c_theta)^nStages - 1);

end

function [y, x, etac, N_target] = operating_point_N_match_drag(altitude, Minf, T_req , To4, gam,R,cp,A2_5,etast,ux,n,alfa1,beta2,rbar)

fn = @(NN) match_drag_subiteration(NN, altitude, Minf, T_req , To4, gam,R,cp,A2_5,etast,ux,n,alfa1,beta2,rbar);

N_target = fzero(fn, 15000);

[~, y, x, etac] = match_drag_subiteration(N_target, altitude, Minf, T_req , To4, gam,R,cp,A2_5,etast,ux,n,alfa1,beta2,rbar);

end

function [fn_value, y, x, etac] = match_drag_subiteration(N, altitude, Minf, T_req , To4, gam,R,cp,A2_5,etast,ux,n,alfa1,beta2,rbar)

U = N/60*2*pi*rbar;
[Ta , aa , Pa , ~] = atmosisa(altitude);

To1 = (1+(gam-1)/2*Minf^2)*Ta;
To2 = To1;
Po1 = Pa*(To1/Ta)^(gam/(gam-1));

% uinf = Minf*sqrt(gam*Rair*Ta);
uinf = Minf * aa;

%Fan properties
bypass = 2.9;	%bypass ratio
FPR = 2.0;

etad = 0.97;
etaf = 0.92;
% etac = 0.87;
etat = 0.91;
etan = 0.98;

Q = 45000e3;

%diffuser
Po2=Pa*(1+etad*(To2/Ta-1))^(gam/(gam-1));

%Fan
Po2_5 = FPR*Po2;
To2_5 = To2*(1 + 1/etaf*((Po2_5/Po2)^((gam-1)/gam) - 1) );
%Calculate compresssor performance based on operating conditions
[y,x, etac] = compmap(gam,R,cp,To2_5,A2_5,etast,ux,U,n,alfa1,beta2);
%y = pressure ratio, P_03 / P_02.5
CPR = y;
mdot_core = x * (0.04*Po2_5*A2_5) / sqrt(To2_5);

% CPR = compfn(cp,To2_5,etast(ux/U),ux,U,n,alfa1,beta2)
Po3=CPR*Po2_5;
To3 = To2_5*(1 + 1/etac*((Po3/Po2_5)^((gam-1)/gam) - 1) );

% Combustor
Po4 = Po3;
f = (To4/To3-1)/(Q/cp/To3-To4/To3);

%turbines: 
%HPT drives compressor
To4_5 = To4 - (To3 - To2_5);
Po4_5 = (1 - 1/etat*(1-To4_5/To4))^(gam/(gam-1))*Po4;

%LPT drives fan
To5 = To4_5 - (1+bypass)*(To2_5 - To2);
Po5 = (1 - 1/etat*(1-To5/To4_5))^(gam/(gam-1))*Po4_5;


%MIX STREAMS!
To6 = ( (1+f)*To5 + bypass*To2_5 ) / ((1+f) + bypass);
Po6 = ( (1+f)*To5 + bypass*To2_5 ) / ((1+f)*To5/Po5 + bypass*To2_5/Po2_5);


%exit condition for mixed stream!
if (Po6 > Pa) 
	Pe = Pa;
	Te = To6 * ( 1 - etan*(1-(Pe/Po6)^((gam-1)/gam)) );
	ue_mix = sqrt( 2*etan * ( gam/(gam-1)*R*To6 * (1-(Pe/Po6)^((gam-1)/gam)) ) );

	%specific thrust: Thrust = mdot_core * TSC
	TSC = (1+f+bypass)*ue_mix - (1+bypass)*uinf;

	T_avail = TSC * mdot_core;
	% mdot_core = Treq/TSC;
	mdot_total = (1+bypass)*mdot_core;

else disp('WHOA!  Mixed-stream pressure lower than ambient, no exit flow... quitting program!')
	T_avail = 0;
	TSFC=1E6;
% 	return;
end

fn_value = T_req - T_avail;
end

function [y x T_avail TSFC] = operating_point(altitude, Minf, T_required , To4, gam,R,cp,A2_5,etast,ux,U,n,alfa1,beta2)

Treq = T_required;

[Ta , aa , Pa , rho_0] = atmosisa(altitude);

To1 = (1+(gam-1)/2*Minf^2)*Ta;
To2 = To1;
Po1 = Pa*(To1/Ta)^(gam/(gam-1));

% uinf = Minf*sqrt(gam*Rair*Ta);
uinf = Minf * aa;

%Fan properties
bypass = 2.9;	%bypass ratio
FPR = 2.0;

etad = 0.97;
etaf = 0.92;
% etac = 0.87;
etat = 0.91;
etan = 0.98;

Q = 45000e3;

%diffuser
Po2=Pa*(1+etad*(To2/Ta-1))^(gam/(gam-1));

%Fan
Po2_5 = FPR*Po2;
To2_5 = To2*(1 + 1/etaf*((Po2_5/Po2)^((gam-1)/gam) - 1) );

%Calculate compresssor performance based on operating conditions
[y,x, etac] = compmap(gam,R,cp,To2_5,A2_5,etast,ux,U,n,alfa1,beta2);
%y = pressure ratio, P_03 / P_02.5
CPR = y;
mdot_core = x * (0.04*Po2_5*A2_5) / sqrt(To2_5);

% CPR = compfn(cp,To2_5,etast(ux/U),ux,U,n,alfa1,beta2)
Po3=CPR*Po2_5;
To3 = To2_5*(1 + 1/etac*((Po3/Po2_5)^((gam-1)/gam) - 1) );

% Combustor
Po4 = Po3;
f = (To4/To3-1)/(Q/cp/To3-To4/To3);

%turbines: 
%HPT drives compressor
To4_5 = To4 - (To3 - To2_5);
Po4_5 = (1 - 1/etat*(1-To4_5/To4))^(gam/(gam-1))*Po4;

%LPT drives fan
To5 = To4_5 - (1+bypass)*(To2_5 - To2);
Po5 = (1 - 1/etat*(1-To5/To4_5))^(gam/(gam-1))*Po4_5;


%MIX STREAMS!
To6 = ( (1+f)*To5 + bypass*To2_5 ) / ((1+f) + bypass);
Po6 = ( (1+f)*To5 + bypass*To2_5 ) / ((1+f)*To5/Po5 + bypass*To2_5/Po2_5);


%exit condition for mixed stream!
if (Po6 > Pa) , Pe = Pa;
else disp('WHOA!  Mixed-stream pressure lower than ambient, no exit flow... quitting program!')
	T_avail = 0;
	TSFC=1E6;
	return;
end

Te = To6 * ( 1 - etan*(1-(Pe/Po6)^((gam-1)/gam)) );
ue_mix = sqrt( 2*etan * ( gam/(gam-1)*R*To6 * (1-(Pe/Po6)^((gam-1)/gam)) ) );

%specific thrust: Thrust = mdot_core * TSC
TSC = (1+f+bypass)*ue_mix - (1+bypass)*uinf;

T_avail = TSC * mdot_core;
% mdot_core = Treq/TSC;
TSFC = f*mdot_core / T_avail * 9.81 * 3600;

mdot_total = (1+bypass)*mdot_core;

end
