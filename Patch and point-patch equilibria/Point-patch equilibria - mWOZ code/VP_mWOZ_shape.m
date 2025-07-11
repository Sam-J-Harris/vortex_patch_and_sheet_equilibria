function [z, zseg, zs, lambda, lambs] = VP_mWOZ_shape(Ra,nsym,csgn,M,b,omega)
% Determines the point vortex - vortex patch configuration
% z = boundary data, zseg = segment data, zs = point vortex locations, 
% lambda = (initial) dimless constant of patch, dimless constant of pt vorts
%
% The initial patch shape is always a circle of radius Ra
% nsym is the symmetry e.g. nsym = 3 => 3 pt vortices
% csgn is the circulation of the vortices (pick either 1 or -1)
% b is the distance of the pt vortex to the origin
%
% Code:
tht = linspace(-pi/2,2*pi/nsym - pi/2, M+1); % theta values for patch segment
zseg = (Ra.*(cos(tht)+1i.*sin(tht))).'; % patch segment "zseg"

thtot = linspace(-pi/2,3*pi/2,(M*nsym)+1); % theta values for full patch
z = (Ra.*(cos(thtot)+1i.*sin(thtot))).'; % full patch "z"

ths = linspace(-pi/2, 3*pi/2, nsym+1); % theta values for satellite pt vortices
zs = (b.*(cos(ths)+1i.*sin(ths))).'; zs = zs(1:end-1); % satellite pt vortices "zs"

A = polyarea(real(z),imag(z)); % dimensional area of initial vortex patch
lambda = A/(2*pi*b^2); % (initial) dimless constant of patch
lambs = csgn/(2*pi*omega*b^2); % dimless constant of pt vorts
end