%% Geometry: 3 equally spaced sheets of unit length emanating from origin
% Method based on Baddoo and Trefethen, Maple Transactions (2021), Appendix B
clc; close all;

M = 200; % number of points on each sheet
s = tanh(linspace(-15,15,M));
omega=1; % fix the angular velocity of the sheet structure
L = 1i;
su = s + 10i*eps; sl = s - 10i*eps;
Z = [-1i*exp(7*1i*pi/6)*(su+1)/2, -1i*L*(su+1)/2, -1i*L*(sl+1)/2,-1i*exp(-1i*pi/6)*(su+1)/2, -1i*exp(-1i*pi/6)*(1+sl)/2, -1i*exp(7*1i*pi/6)*(1+sl)/2].';
z1 = -1i*exp(7*pi*1i/6); z2 = 0; z3 = -1i*L; z4 = -1i*exp(-1i*pi/6); % vortex locations (only 3 vortices here so z2=0).

%% Joukowsky variable with branch cut [a,b]
w = @(z,a,b) z - (a+b)/2 + exp((log(exp(-1i*angle(a-b))*(z-a)) ...
+ log(exp(-1i*angle(a-b))*(z-b)))/2 ...
+ 1i*angle(a-b));

%% Matrix problem
nL = 0; nP =20; s0 = nP/2;
lc = @(z,p,a,b) 1./(s0-log(1-w(p,a,b)./w(z,a,b)));
[H1a,P1 ] = VAorthog(1./w(Z,z1,z2),nL); % segment 1
[H1b,L1b] = VAorthog(lc(Z,z1,z1,z2),nP);
[H1c,L1c] = VAorthog(lc(Z,z2,z1,z2),nP);
[H2a,P2 ] = VAorthog(1./w(Z,z2,z3),nL); % segment 2
[H2b,L2b] = VAorthog(lc(Z,z2,z2,z3),nP);
[H2c,L2c] = VAorthog(lc(Z,z3,z2,z3),nP);
[H3a,P3 ] = VAorthog(1./w(Z,z2,z4),nL); % segment 3
[H3b,L3b] = VAorthog(lc(Z,z2,z2,z4),nP);
[H3c,L3c] = VAorthog(lc(Z,z4,z2,z4),nP);
A = [P1 L1b L1c P2 L2b L2c P3 L3b L3c]; 
A = [real(A) -imag(A)]; % matrix of basis vectors
source = @(z)  -(omega/2)*abs(z).^2+(.5797635)*log(w(z,z1,z2)/2); % choose coeff .57972 by trial and error to obtain stagnation points at tips
F = real(source(Z));
cc = A\F;

%% Function handles for solution
c = reshape(cc,[],2)*[1; 1i];
g = @(z) reshape([VAeval(1./w(z(:),z1,z2),H1a),...
VAeval(lc(z(:),z1,z1,z2),H1b)...
VAeval(lc(z(:),z2,z1,z2),H1c)...
VAeval(1./w(z(:),z2,z3),H2a),...
VAeval(lc(z(:),z2,z2,z3),H2b)...
VAeval(lc(z(:),z3,z2,z3),H2c)...
VAeval(1./w(z(:),z2,z4),H3a)...
VAeval(lc(z(:),z2,z2,z4),H3b)...
VAeval(lc(z(:),z4,z2,z4),H3c)]*c,size(z));
Gcmplx = @(z) source(z) - g(z);
G = @(z) real(Gcmplx(z)); H = @(z) imag(Gcmplx(z));

%% Contour plot
figure(1) % Slit and streamlines
LW = 'LineWidth'; PO = 'position';
% axes(PO,[.03 .36 .57 .57])
x = linspace(-2,2,500); y = linspace(-2,2,500);
[xx,yy] = meshgrid(x,y); zz = xx+1i*yy;
HH = H(zz); HH(G(zz)>1) = NaN;
HH(real(zz)<-1 & abs(imag(zz))<0.1) = NaN;
levelsH = linspace(min(min(HH)),max(max(HH)),14); levelsH([1,end]) = [];
levels = -.5:.05:.5; contour(x,y,G(zz),levels,'k'), clim([0 1.1]), hold on
plot([z1 z2 z3 z2 z4],'b',LW,1.5), axis([-1.7 1.7 -1.7 1.7]), axis square, axis off, hold on
hold off

figure(2) % Circulation density along one of the sheets
h=0.00001;
X=h:h:1-h;
phi_sheet=imag(Gcmplx(X));
sheetvort=-2*diff(phi_sheet)/h;
vel_tip=sheetvort(end); % iterate 'coeff' above to make this value close to zero
fprintf("Velocity at sheet tip = " + num2str(vel_tip)+ ".\n")
plot(X(:,1:length(sheetvort)),sheetvort)
xlabel('$x$','interpreter','latex','FontSize',18)
ylabel('$\rho(x)$','interpreter','latex','FontSize',18)
set(gca,'xtick',0:.2:1,'ytick',0:.2:1.4,'fontsize',14)
total_circ=trapz(X(:,1:length(sheetvort)),sheetvort);

hold off