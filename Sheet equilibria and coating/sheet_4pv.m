function sheet_4pv % Stream function exterior to a vortex sheet in rotating frame with angular speed omega
%
% We seek the stream function u(z) that is zero on the interval c + r[-1,1] and harmonic outside
% this slit, including at z=infty, except with u(z)\simlog|z-v_c| as z->v_c are the 4 vortex locations. u is expanded as
% u(z) = .5*omega*|z|^2 + (sheet_circulation)*log|w| + a(1)
%           + SUM_{k=1}ˆN a(2k)*real(wˆ(-k)) + a(2k+1)*imag(wˆ(-k))
%
% where w(z) is a Joukowski map of exterior(slit) to exterior(unit disk).
% Given total sheet circulation ('sheet'_circulation), 'fsolve' finds omega, 2*r and the location of the vortices on the imaginary axis and plots (i) streamlines in rotating frame and 
% (ii) circulation density along sheet
%
%% Main Code
clc; close all;
format long

N = 80; npts = 3*N; % no. expansion terms and sample pts
circ = (1+1e-12)*exp(2i*pi*(1:npts)'/npts);
sheet_circulation=-1; % units of 2pi
vortloc=1; % vortex location on real axis
vortloc2_guess=.9005*1i; % vortex location on imaginary axis (initial guess)
gamma=-1; % vortex strength (units of 2pi)
c=0; % centre of slit
r_guess=0.75; % half length of the slit (initial guess)
omega_guess=2.0; % angular speed (initial guess)
guess=[omega_guess, r_guess, vortloc2_guess];
options = optimoptions('fsolve','Display','off');
sol=fsolve(@eqns,guess, options);

omega=sol(1); r=sol(2); vortloc2=sol(3); % angular velocity, half length of slit, vortex location (on imaginary axis)
startslit=c-r; endslit=c+r; % start and end points of slit
a=lsq(sol); % complex coefficients from least-squares method

fprintf("Solutions: omega = "+num2str(omega)+", r = "+num2str(r)+", vortex location (on imaginary axis) = "+num2str(vortloc2)+".\n")

fprintf("Real axis vortex velocity = " + num2str(dzdt_vortloc(sol,vortloc,a))+".\n")
fprintf("Imaginary axis vortex velocity = " + num2str(dzdt_vortloc2(sol,vortloc2,a))+".\n")

%% Figures
figure(5); % Slit and streamlines
x = linspace(-3.,3,1000); y = linspace(-3,3,1000); [xx,yy] = meshgrid(x,y);
zz = xx+1i*yy; uu = slit1fun(zz); % uu = stream function
z = .00001*1i+c + r*[-1 1];
plot(z,'b','linewidth',1.6), hold on
levels = -2:.25:1.25; contour(xx,yy,uu,levels,'k'), axis equal, axis([-2.25 2.25 -2.25 2.25]), axis off
set(gca,'xtick',-3:1:3,'ytick',-3:1:3,'fontsize',8)
hold off 

figure(6); % Circulation density along the sheet
eps=0*.000001;
xxx=linspace(startslit+eps,endslit-eps,2000)+0*.0000001*1i+1e-12*1i; % coordinates along the slit (with distance eps from start/end of slit)
compvelplot=dzdt_sheet(xxx,a); % grad(psi)
compvelplot(1)=0;
compvelplot(2000)=0;
plot(real(xxx),2*imag(compvelplot)); % circulation density vs x-position along the slit
axis([-0.8 0.8 0 6])
xlabel('$x$','interpreter','latex','FontSize',18)
ylabel('$\rho(x)$','interpreter','latex','FontSize',18)
set(gca,'xtick',-.8:.4:.8,'ytick',0:2:6,'fontsize',14)
hold on

% % O'Neil exact solution (uncomment if desired)
% yy=linspace(-1,1,1000);
% TT=5*gamma*2*pi;
% circ2=gamma*2*pi;
% xv1=1;
% xv2=-1;
% yv1=vortloc2;
% yv2=-vortloc2;
% om=-omega;
%
% plot(yy,-2*real(sqrt(-omega^2*yy.^2+om*TT/pi-(circ2^2/(4*pi^2))*(1./(yy-xv1).^2+1./(yy-xv2).^2+1./(yy-yv1).^2+1./(yy-yv2).^2)+(om*circ2/pi)*((yv1-yv2)./(yy-yv1)+(yv2-yv1)./(yy-yv2)) )))
hold off 

%% Additional functions
function F=eqns(x) % equation function for fsolve
a=lsq(x);
F(1)=sheetprop(x,a)-2*pi*sheet_circulation; % angular speed
F(2)=abs(dzdt_vortloc(x,vortloc,a)); % half length of slit
F(3)=abs(dzdt_vortloc2(x,x(3),a)); % vortex location 2
end

function w = wz(x,z)
zc = (z-c)/x(2); sgn = real(zc)>0|(real(zc)==0&imag(zc)>0); sgn = 2*sgn - 1;
w = zc + sgn.*sqrt(zc.^2-1);
end

function coeffs=lsq(x) % least-squares function
z = c + x(2)*(circ+1./circ)/2; % sample points - note x(2) from sol is r = half length of the slit
w = wz(x,z); 
rhs = -x(1)*0.5*abs(z).^2-gamma*log(abs(z-vortloc))-gamma*log(abs(z+vortloc))-gamma*log(abs(z-1i*imag(x(3))))-gamma*log(abs(z-1i*imag(x(3)))); % right-hand side
A = ones(npts,2*N+1);
for k = 1:N % set up least-squares matrix
A(:,2*k ) = real(w.^-k);
A(:,2*k+1) = imag(w.^-k);
end
coeffs = A\rhs; % complex coefficient vector
end 

function G=sheetprop(x,a)
startslit=c-x(2);
endslit=c+x(2);
xxx=linspace(startslit,endslit,1000)+.0000001*1i;
compvel=dzdt(x,xxx,a);
G=2*(trapz(real(xxx),-imag(compvel)));% vorticity of sheet is the velocity jump across the sheet
end

function GG=sheetprop_cent(x,a)
startslit=c-x(2);
endslit=c+x(2);
xxx=linspace(startslit,endslit,1000)+.0000001*1i;
compvel=dzdt(x,xxx,a);
GG=2*(trapz(real(xxx),-real(xxx).*imag(compvel)));% centre of vorticity of sheet 
end

function g = dzdt(x,z,a)
zc = (z-c)/x(2); sgn = real(zc)>0|(real(zc)==0&imag(zc)>0); sgn = 2*sgn - 1;
w = zc+sgn.*sqrt(zc.^2-1); dwdz = (1+sgn.*zc./sqrt(zc.^2-1))/x(2);
g = x(1)*conj(z)+gamma./(z-vortloc)+gamma./(z+vortloc)+gamma./(z-1i*imag(x(3)))+gamma./(z+1i*imag(x(3)))+sheet_circulation*dwdz./w;
for k = 1:N, g = g - k*(a(2*k)-1i*a(2*k+1))*dwdz./w.^(k+1); end
end

function g = dzdt_vortloc(x,z,a) % velocity of real axis vorticities (should be close to 0 in equilibrium)
zc = (z-c)/x(2); sgn = real(zc)>0|(real(zc)==0&imag(zc)>0); sgn = 2*sgn - 1;
w = zc+sgn.*sqrt(zc.^2-1); dwdz = (1+sgn.*zc./sqrt(zc.^2-1))/x(2);
g = x(1)*conj(z)+gamma./(z+vortloc)+gamma./(z-1i*imag(x(3)))+gamma./(z+1i*imag(x(3)))+sheet_circulation*dwdz./w;
for k = 1:N, g = g - k*(a(2*k)-1i*a(2*k+1))*dwdz./w.^(k+1); end
end

function g = dzdt_vortloc2(x,z,a) % velocity of imaginary axis vorticities (should be close to 0 in equilibrium)
zc = (z-c)/x(2); sgn = real(zc)>0|(real(zc)==0&imag(zc)>0); sgn = 2*sgn - 1;
w = zc+sgn.*sqrt(zc.^2-1); dwdz = (1+sgn.*zc./sqrt(zc.^2-1))/x(2);
g = x(1)*conj(z)+gamma./(z+vortloc)+gamma./(z-vortloc)+gamma./(z+x(3))+sheet_circulation*dwdz./w;
for k = 1:N, g = g - k*(a(2*k)-1i*a(2*k+1))*dwdz./w.^(k+1); end
end

% Functions for plotting are below
function u = slit1fun(z)
a=lsq(sol);
w = wz(sol,z); u = 0.5*omega*abs(z).^2+gamma*log(abs(z-vortloc))+gamma*log(abs(z+vortloc))+gamma*log(abs(z-vortloc2))+gamma*log(abs(z+vortloc2))+sheet_circulation*log(abs(w)) + a(1);
for k = 1:N, u = u + a(2*k)*real(w.^(-k)) + a(2*k+1)*imag(w.^(-k)); end
end

function sheetvel = dzdt_sheet(z,a) % grad(psi) = complex velocity on sheet
w = wz(sol,z); 
zc = (z-c)/r; sgn = real(zc)>0|(real(zc)==0&imag(zc)>0); sgn = 2*sgn - 1;
w = zc+sgn.*sqrt(zc.^2-1); dwdz = (1+sgn.*zc./sqrt(zc.^2-1))/r;
sheetvel = omega*1*conj(z)+gamma./(z-vortloc)+gamma./(z+vortloc)+gamma./(z-vortloc2)+gamma./(z+vortloc2) +sheet_circulation*dwdz./w;
for k = 1:N, sheetvel = sheetvel - k*(a(2*k)-1i*a(2*k+1))*dwdz./w.^(k+1); end
end

end

