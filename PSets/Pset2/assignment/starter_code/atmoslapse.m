function [T, a, P, rho] = atmoslapse(h, g, gamma, R, L, hts, htp, rho0, P0, T0 )
%  ATMOSLAPSE Use Lapse Rate Atmosphere Model.
%   [T, A, P, RHO] = ATMOSLAPSE(H, G, GAMMA, R, L, HTS, HTP, RHO0, P0, T0) 
%   implements the mathematical representation of the lapse rate atmospheric 
%   equations for ambient temperature, pressure, density, and speed of sound 
%   for the input geopotential altitude.  This atmospheric model is customizable
%   by specifying the atmospheric properties in the function input.
%
%   Inputs required by ATMOSLAPSE are:
%   H      :a numeric array of M geopotential height in meters. 
%   G      :a scalar of acceleration due to gravity in meters per second squared.
%   GAMMA  :a scalar of specific heat ratio.
%   R      :a scalar of characteristic gas constant joules per kilogram-kelvin.
%   L      :a scalar of lapse rate in kelvin per meter.
%   HTS    :a scalar of height of troposphere in meters.
%   HTP    :a scalar of height of tropopause in meters.
%   RHO0   :a scalar of air density at mean sea level in kilograms per meter cubed.
%   P0     :a scalar of static pressure at mean sea level in pascal.
%   T0     :a scalar of absolute temperature at mean sea level in kelvin.
%
%   Outputs calculated for the lapse rate atmosphere are: 
%   T      :a numeric array of M temperature in kelvin.
%   a      :a numeric array of M speed of sound in meters per second.
%   P      :a numeric array of M air pressure in pascal.
%   rho    :a numeric array of M air density in kilograms per meter cubed.
%
%   Limitation: 
%
%   Below the geopotential altitude of 0 km and above the geopotential
%   altitude of the tropopause, temperature and pressure values are held.
%   Density and speed of sound are calculated using a perfect gas
%   relationship.
%
%   Example:
%
%   Calculate the atmosphere at 1000 meters with the International Standard 
%   Atmosphere input values:
%      [T, a, P, rho] = atmoslapse(1000, 9.80665, 1.4, 287.0531, 0.0065, ...
%          11000, 20000, 1.225, 101325, 288.15 );
%
%   See also ATMOSCIRA, ATMOSCOESA, ATMOSISA, ATMOSNONSTD.

%   Copyright 2000-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $  $Date: 2010/09/28 03:13:37 $

%   Reference:  U.S. Standard Atmosphere, 1976, U.S. Government Printing 
%   Office, Washington, D.C.

error(nargchk(10,10,nargin,'struct'));

if ~(isscalar(g)&&isscalar(gamma)&&isscalar(R)&&isscalar(L)&&isscalar(hts) ...
     &&isscalar(htp)&&isscalar(rho0)&&isscalar(P0)&&isscalar(T0))
    error(message('aero:atmoslapse:nonScalar'));
end

if ~isnumeric( h )
    error(message('aero:atmoslapse:notNumeric'));
end

for i = length(h): -1 :1
    if ( h(i) > htp )
        h(i) = htp;
    end

    if ( h(i) <  0 )
        h(i) = 0;
    end

    if ( h(i) > hts )
        T(i) = T0 - L*hts; 
        expon(i) = exp(g/(R*T(i))*(hts - h(i))); 
    else
        T(i) = T0 - L*h(i); 
        expon(i) = 1.0;  
    end
end

a = sqrt(T*gamma*R);

theta = T/T0;

P = P0*theta.^(g/(L*R)).*expon;
rho = rho0*theta.^((g/(L*R))-1.0).*expon;
