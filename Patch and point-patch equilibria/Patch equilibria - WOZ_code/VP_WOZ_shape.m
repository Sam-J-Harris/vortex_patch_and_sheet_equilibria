function [z, zs, Ra, A] = VP_WOZ_shape(nsym,shswt,M,Pts) % shape data - rotating V-state
% z = boundary data, zs = segment data, Ra = radius distance of start segment point, A = area of patch
Ratab = [[1.05 1.22128 1.39256 1.56384 1.73512]; [1.05 1.14709 1.24418 1.34127 1.43836]; [1.05 1.11598 1.18174 1.24761 1.31348]; [1.05 1.09855 1.14711 1.19566 1.24421]]; % table of Ra values: those used in Wu et al. 1984 - see their Table III 
Ra = Ratab(max(nsym-2,1),shswt); % for given nfold symmetry (= 3, 4, 5 or 6), picks Ra from table based on shswt

tht = linspace(-pi/2,pi/nsym - pi/2, M+1); % theta values for segment
R = Ra + (1-Ra).*(tht+pi/2)./(pi/nsym); % radial values for segment
zs = (R.*(cos(tht)+1i.*sin(tht))).'; % segment data

thtot = linspace(-pi/2,3*pi/2,Pts); % full theta values
Rseg = [R flip(R(2:end-1))]; % gives R values on one "fold"
Rtot = Rseg; 
for k=1:nsym-1
    Rtot = [Rtot Rseg]; % add segments
end
Rtot = [Rtot Rtot(1)]; % full R values
z = (Rtot.*(cos(thtot)+1i.*sin(thtot))).'; % full boundary data
A = polyarea(real(z),imag(z)); % dimensional area
end