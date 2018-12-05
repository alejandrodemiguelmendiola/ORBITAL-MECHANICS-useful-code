function dy = ode_2bodyproblem(~,y,mu)
%ODE defining the orbital dynamics.
%
%PROTOTYPE
%   dy = ode_2bodyproblem(t,y,mu)  
%
%INPUT
%   t[1]    Time 
%   y[1x6]  State vector [rx ry rz vx vy vz]
%   mu[1]   Gravitational parameter (constant) 
%
%OUTPUT
%   dy[6x1] Derivative of the state

%Module of the radiovector.
r = sqrt(y(1)^2 + y(2)^2 + y(3)^2);

%Set the derivatives of the state.
dy = [
    y(4);
    y(5);
    y(6);
    -mu/(r^3) * y(1);
    -mu/(r^3) * y(2);
    -mu/(r^3) * y(3)
    ]; 
end
