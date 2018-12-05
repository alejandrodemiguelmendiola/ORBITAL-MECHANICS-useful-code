function [rp,delta,Deltavp] = poweredGA(vm_inf,vp_inf,mu_planet,radius_planet,h_atm_planet)
%
% PROTOTYPE [rp,delta,Deltavp] = poweredGA(vm_inf,vp_inf,mu_planet,radius_planet,h_atm_planet)
% 
% INPUT:
%   vm_inf[3x1]      incoming velocity of hyperbola
%   vp_inf[3x1]      outcoming velocity of hyperbola
%   mu_planet
%   radius_planet
%   h_atm_planet
% 
% OUTPUT:
%   rp[1x1]        perigee of the hyperbola 
%   delta[1x1]     turning angle [rad]
%   Deltavp[1x1]   required Deltav for the manouvre

%% Compute the turning angle of the GA
delta = acos(dot(vp_inf,vm_inf)/(norm(vm_inf)*norm(vp_inf))); %rad

%% Solve the non-linear problem to compute rp
initguess = radius_planet;
opts = optimoptions('fsolve','OptimalityTolerance',1e-13);
rp = fsolve(@(x) root(x,norm(vm_inf),norm(vp_inf),delta,mu_planet),initguess,opts);

%% Check rp validity
validity = rp > radius_planet + h_atm_planet;

%% Compute other useful variables (inc/outc eccentricity and turning angle)
em = (1+rp*norm(vm_inf)^2/mu_planet);
ep = (1+rp*norm(vp_inf)^2/mu_planet);
deltam = 2*asin(1/em);
deltap = 2*asin(1/ep);

%% Compute velocities at pericenter and Deltavp
vmp_inf = sqrt(norm(vm_inf)^2 + 2*mu_planet/rp); %eq.8.58 from Curtis, v at perigee
vpp_inf = sqrt(norm(vp_inf)^2 + 2*mu_planet/rp);
Deltavp = vpp_inf - vmp_inf;
h_ga = rp - radius_planet;

%% Defining the FSOLVE function
function F = root(x,vm_inf,vp_inf,delta,mu_planet)

    rp = x(1);
    em = (1+rp*(vm_inf)^2/mu_planet);
    ep = (1+rp*(vp_inf)^2/mu_planet);
    deltam = 2*asin(1/em);
    deltap = 2*asin(1/ep);
    
    F(1) = delta - deltam/2-deltap/2;
end


end

