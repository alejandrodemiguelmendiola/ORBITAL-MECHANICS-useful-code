function [T] = getT(r,v,mu)
%   Calculates the period of an orbit from r,v and mu
% PROTOTYPE: [T] = getT(r,v,mu)
% 
% INPUT:
%       r[1x3]          Position vector of the orbit [km]
%       v[1x3]          Velocity vector of the orbit [km/s]
%       mu[1x1]         Gravitational parameter [km^3/s^2]
% 
% OUTPUT:
%       T[1x1]          Orbit Period [s]

E = norm(v)^2/2-mu/norm(r); % Specific Energy
a = -mu/2/E; % Semi major axis
n = sqrt(mu/a^3); % Mean velocity
T = abs(2*pi/n); % Period
end

