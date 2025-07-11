% Solve coating problem surrounding finite substrate
clear; clc; close all; tic; % close all, start timer
options = optimoptions('fsolve','Display','iter');
warn = warning('off','MATLAB:rankDeficientMatrix');  % matrices are ill-conditioned

M=2000; % Number of points on coating layer and substrate 
theta=linspace(0,2*pi*(M-1)/M,M); zeta=exp(1i*theta);
A=3; % set initial radius of circular coating boundary
rr=A*abs(zeta);
zz0=rr.*zeta; % z coordinates of initial circular coating boundary
axlim=3;
plot(real(zz0), imag(zz0),'r') % plot inital coating shape guess in red
axis([-axlim axlim -axlim axlim]), axis square
hold on

% Substrate geometry
c2=0;
b=1;% half-height of block
hhh=linspace(-1,1,M/4); % inner Jordan curve
zsubs = 0.5*[hhh-b*1i 1i*(b*hhh-1*1i) hhh+b*1i -1i*(b*hhh-1*1i)].'; % z coordinates of substrate boundary
fprintf("Area enclosed by substrate: "+num2str(polyarea(real(zsubs),imag(zsubs)))+"\n") % area enclosed by the substrate

% Boundary condition on substrate (positive number)
boundvalue=.02;

% Perform iterations with relaxation parameter to obtain boundary with zero gradient
for k=1:25
zz=rr.*zeta;
% Obtain first and second derivative of the harmonic part of solution
[first,second]=gradbound(zz,M,zsubs,boundvalue);
% Construct variation of radii: delta r=-u_r/u_rr
del=0*1./abs(zz)+0.5*abs(zz)+real(zz.*first)./abs(zz);
sec_deriv=-0*1./abs(zz).^2+0.5+real(zz.^2.*second)./(abs(zz).^2);
% Stopping condition
tolstep=rms(del)/rms(rr);
fprintf("Tolstep at k="+num2str(k)+": " + num2str(tolstep)+"\n")
if tolstep < 10^(-4) %stopping criteria
    fprintf("Total iterations: " +num2str(k)+"\n");
break
end

relax=.8;% Relaxation parameter
rr=rr-relax*del./sec_deriv; % change in r based on relaxation param
zz=rr.*zeta;
plot(real(zz), imag(zz),'--')
end
hold off
figure
% Plot final shape (black solid line)
zz=rr.*zeta;
plot(real(zz), imag(zz),'k')
axis([-1.2 1.2 -1.2 1.2])
axis square
axis off
hold on

% Plot substrate
plot(real(zsubs),imag(zsubs),'b','LineWidth',2.0)
hold off

%------------------------------------------------------------------------------------------------------

function [F,G]=gradbound(zz,N,z2,subsbc)
boundaryconst=subsbc; % value of boundary constant
%% Boundaries
c1=0;c2=0;
%M=800; %number of points

z1 = zz.';
Z = [z1; z2]; % both coating and substrate boundaries (for generating poles)

% Construct BCs for psi=phi+r^2/4, where phi is harmonic
hh=@(z) -(abs(z).^2)/4-0*log(abs(z));
H=hh(Z);
for j=1:N
H(N+j)=boundaryconst+H(N+j);
end

%% AAA Algorithm and Poles
pol_in = []; pol_in2 = [];
tol=1e-8;

% Global AAA
[~,polk] = aaa(H,Z,'tol',tol,'cleanup',0); titlegl='Global'; % global AAA
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); % determines if point is in a polygon
polk_in = polk(~inpolygonc(polk,z1)|inpolygonc(polk,z2)); pol_in = [pol_in; polk_in];
polk_in2 = polk(inpolygonc(polk,z2)|~inpolygonc(polk,z1)); pol_in2 = [pol_in2; polk_in2];

%% Least Squares
pol1 = pol_in.'; pol2 = pol_in2.'; 
plot(real(pol1),imag(pol1),'.k')
plot(real(pol2),imag(pol2),'.r')
d1 = min(abs((z1)-pol1),[],1); d2 = min(abs((z2)-pol2),[],1);
K=20;

P1=(Z-c1).^(0:K);
Q1 = d1./((Z)-pol1);

P2 = (1./(Z-c2)).^(0:K);
Q2 = d2./((Z)-pol2);
A = [0*Z log(abs((Z-c2))) real(P1) real(Q1) real(P2) real(Q2) 0*Z 0*Z -imag(P1) -imag(Q1)  -imag(P2) -imag(Q2)];
c = reshape(A\H,[],2)*[1; 1i];

% Funtion f and its first (fd) and second (fdd) derivatives
f = @(z) reshape([0*z(:) log((z(:)-c2)) (z(:)-c1).^(0:K) d1./(z(:)-pol1) 1./(z(:)-c2).^(0:K) d2./(z(:)-pol2)]*c,size(z));
fd= @(z) reshape([0*z(:) 1./(z(:)-c2) (0:K).*(z(:)-c1).^((0:K)-1) -d1./((z(:)-pol1).^2) -(0:K)./(z(:)-c2).^(1:K+1) -d2./((z(:)-pol2).^2)]*c,size(z));
fdd= @(z) reshape([0*z(:) -1./(z(:)-c2).^2 (0:K).*((0:K)-1).*(z(:)-c1).^((0:K)-2) 2*d1./((z(:)-pol1).^3) (0:K).*(1:K+1)./(z(:)-c2).^(2:K+2) 2*d2./((z(:)-pol2).^3)]*c,size(z));

u = @(z) -hh(z)+real(f(z));%Laplace solution is phi=real(f); Poisson solution u=phi+(r^2)/4
gradu = @(z) fd(z); % note gradu is complex-valued

% zero gradient condition at outer boundary
for j=1:N
F(j)=gradu(zz(j)); % first derivative 
G(j)=fdd(zz(j)); % second derivative
end

end
