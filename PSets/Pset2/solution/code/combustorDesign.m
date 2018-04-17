function L = combustorDesign(UInf, altitude, engine, fluid)
%This outputs the zonal lengths of a nominal RQL combustor.
    %get the properties
    [T0,a0,p0,~] = atmosisa(altitude);
    %given
    T04=engine.turbineInletTemperature; %[K]
    pic = engine.overallCompressionRatio; %[]
    cp=fluid.cp; % [J/kg]
    R=fluid.R; % [J/kg.K]
    M0 = UInf/a0;
    LHV = fluid.LHV;
    %compute quantities
    %2.1.1
    g=cp/(cp-R);
    taur = 1+(g-1)/2*M0^2;
    tauc = pic^((g-1)/g);
    T03 = T0*taur*tauc; %[K]
    f = (T04-T03)/(LHV/cp-T04);
    %2.1.2
    phi = f/fluid.fst;
    %2.1.3
    pir = taur^(g/(g-1));
    pic = tauc^(g/(g-1));
    p03 = p0*pir*pic;
    r03 = p03/(R*T03);
    %2.2.a
    Ro = engine.combustor.Ro; %given [m]
    Ri = engine.combustor.Ri; %given [m]
    Rc = (Ro-Ri)/2;
    theta = asin(Rc/(Ri+Rc));
    nComb = floor(2*pi/(2*theta));
    %2.2.b
    Ac = pi*Rc^2;
    %2.3.b.i
    mDotA = engine.mDotA; % 
    betaC = engine.combustor.betaC; 
    nBetaC = length(betaC);
    betaMix = 100e-6;
    Ua1 = mDotA./((1+betaC).*(1-betaMix).*r03.*Ac.*nComb); %[m/s]
    %2.3.b.ii
    fMix = f.*(1+betaC);
    phiMix = phi.*(1+betaC);
    %2.3.b.iii
    D0 = engine.combustor.D0;
    beta = engine.combustor.beta;
    tauEvap = D0.^2.0./beta;
    %2.3.b.iv
    Ue = engine.combustor.Ue;
    L.LEvap = Ue.*tauEvap.*ones(size(betaC));
    %2.3.b.v
    TFlame = zeros(1,nBetaC);
    TEvap = T03; %isothermal
    gas = fluid.gas;
    for iBetaC=1:nBetaC
        X = zeros(nSpecies(gas),1);
        X(speciesIndex(gas,'nc12h26'))=phiMix(iBetaC)/2.0; %only half is combusted
        X(speciesIndex(gas,'o2'))=37/2;
        X(speciesIndex(gas,'n2'))=37/2*3.76;
        set(gas,'T',T03,'P',p03,'X',X);
        equilibrate(gas,'HP');
        TFlame(iBetaC) = temperature(gas);
    end
    UFlame = Ue.*(TFlame./TEvap);
    %2.3.e.i
    dDil = engine.combustor.dDil; 
    nDil = engine.combustor.nDil; 
    UJet = 4.*mDotA.*betaC./((1+betaC).*nComb.*nDil.*pi.*r03.*dDil.^2);
    %2.3.d.ii
    cf =betaC./((1+betaC).*(1+f));
    %2.3.e.iii
    CDil = engine.combustor.CDil; 
    rFlame = p03./(R*TFlame);
    xi = sqrt(r03.*UJet./(rFlame.*UFlame));
    L.LDil = (cf/CDil).^(-3/2).*xi.*dDil;
    %2.3.e.iv
    rDil=r03.*(((1+betaC).*(1+f))./((1+(1+betaC).*f).*TFlame./T03+betaC));
    UDil = mDotA.*(1+f)./(nComb.*rDil.*Ac);
    TDil = TFlame.*(rFlame./rDil);
    %2.3.f.1 
    Ta = 2100; %given [K]
    A = 240; %given [Pa.s]
    n = 1.0; %given
    m = 0.2; %given
    TSrz = TDil;
    phiSrz = phi*(1+betaC)./(1+2.0*betaC);
    tauSrz = A./(p03.^n.*phiSrz.^m).*exp(Ta./TSrz);
    L.LSrz = UDil.*tauSrz;

    %total length
    L.LTotal = L.LEvap+L.LDil+L.LSrz;
end

