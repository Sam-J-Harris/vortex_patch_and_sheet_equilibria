function [psia,ubdry,uvort] = VP_mWOZ_aaa(Z,zseg,zs,lambda,lambs)
% AAA-LS numerical method
%
% psia = total stream function - only needed when plotting s'lines
% ubdry = u+iv = velocity on each boundary point of the patch
% uvort = u+iv = velocity at a point vortex*
%
% *by symmetry, all point vortices have the same velocity magnitude
%
% Code:
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); % determines if point is in polygon (BOUNDARY IS INCLUDED IN POLYGON)

% Normal vector
dZ = circshift(Z(1:end-1),1) - circshift(Z(1:end-1),-1); v = 1i.*(dZ)./abs(dZ); 
nvec = v./abs(v); nvec(end+1) = nvec(1); % unit normal vector

% Far field condition into Dirichlet and Neumann BCs - see Sec. 3
hd = @(z) lambda.*log(abs(z)) - 0.25.*(abs(z).^2); % Dirichlet BC 
hn = @(z) lambda.*real(nvec./z) - 0.5.*(real(nvec.*conj(z))); % Neumann BC
Hd = hd(Z); Hn = hn(Z); bigH = [Hd; Hn]; % apply BCs to each point on the curve; put both into big vector "bigH"

% Global AAA poles %%%%%%%% 
rtoly = 1e-8; ztoly = 1e-2; %rtoly = 0; ztoly = 0; - use these from no pole control
[~,poltot1,restot1] = aaa(Hd,Z,'cleanup',1,'toly',rtoly); % poles for Dirichlet BC
[~,poltot2,restot2] = aaa(Hn,Z,'cleanup',1,'toly',rtoly); % poles for Neumann BC

poltot = [poltot1; poltot2].'; restot = [restot1; restot2].'; % all poles and residues

% Pole control - manual elimination of Froissart doublets  %%%%%%%%
poltot = poltot(abs(restot)>rtoly); % poles of small residue < rtoly eliminated
D = min(abs(Z-poltot),[],1); D(D<ztoly*max(abs(Z))) = 0; D(D~=0) =1; % finds if min dist. b/ween poles and bdry is below toly
poltot=poltot.*D; poltot(abs(poltot)==0)=[]; % poles close to the boundary < ztoly eliminated

polint = (poltot(~inpolygonc(poltot,Z))); % exterior poles (used in interior pb)
polext = (poltot(inpolygonc(poltot,Z))); % interior poles (used in exterior pb)

% % Plot of poles (OPTIONAL) %%%%%%%%
% figure(3) % uncomment if you want to see figure of the poles
% zplot = Z(1:end); xp = real(zplot); yp = imag(zplot);
% plot(xp,yp,'k'), hold on, plot(polint,'.r'); plot(polext,'.b'), hold off, 
% axis square, daspect([1 1 1]); axis([-10 10 -10 10]);

% Interior and exterior poles with min distance to bdry
if size(polint,1)==0 && size(polint,2)==0 % interior pb poles - if no poles were generated
    polint=0; dint=0;
else
    dint = min(abs((Z)-polint),[],1); % dint = distance between interior pb pole and closest point on polygon
end

if size(polext,1)==0 && size(polext,2)==0 % exterior pb poles - if no poles were generated
    polext=0; dext=0; 
else
    dext = min(abs((Z)-polext),[],1); % dext = distance between exterior pb pole and closest point on polygon
end

% Least-Squares %%%%%%%% 
N = 20; N2int = size(polint,2); % N = polynomial series truncation, N2int = number of interior poles (used later when separating interior and exterior coefficients)
P1d = Z.^(0:N); % interior Dirichlet Runge bases
P2d = (1./Z).^(1:N); % exterior Dirichlet Runge bases (Dirichlet BC, polynomial part)

P1n = (0:N).*Z.^(-1:N-1).*nvec; % interior Neumann Runge bases
P2n = -(1:N).*Z.^(-(2:N+1)).*nvec; % exterior Neumann Runge bases (Neumann BC, polynomial part)

Q1d = dint./((Z-polint)); % interior Dirichlet Newman bases
Q2d = dext./((Z-polext)); % exterior Dirichlet Newman bases (Dirichlet BC, poles part)

Q1n = -dint./((Z-polint).^2).*nvec; % interior Neumann Newman bases
Q2n = -dext./((Z-polext).^2).*nvec; % exterior Neumann Newman bases (Neumann BC, poles part)
      
A1d = [real(P1d) real(Q1d) -imag(P1d) -imag(Q1d)]; % matrix of interior Dicichlet bases
A2d = [real(P2d) real(Q2d) -imag(P2d) -imag(Q2d)]; % matrix of exterior Dirichlet bases

A1n = [real(P1n) real(Q1n) -imag(P1n) -imag(Q1n)]; % matrix of interior Neumann bases
A2n = [real(P2n) real(Q2n) -imag(P2n) -imag(Q2n)]; % matrix of exterior Neumann bases

bigA = [A1d A2d; A1n A2n]; testc = bigA\bigH; % combined matrix of basis vectors "bigA"; vector of all coefficients

c1 = testc(1:(2*(N+1)+2*N2int),:); c2 = -testc((2*(N+1)+2*N2int)+1:end,:); % c1 = coefficients of interior terms: polynomial and the N2int poles. c2 = coefficients of exterior terms
c1 = reshape(c1,[],2)*[1;1i]; c2 = reshape(c2,[],2)*[1;1i]; % reshape to be complex coefficients

% Stream function psi
F1 = @(z) reshape([z(:).^(0:N) dint./(z(:)-polint)]*c1, size(z)); %analytic function F1(z)
F2 = @(z) reshape([(1./z(:)).^(1:N) dext./(z(:)-polext)]*c2, size(z)); %analytic function F2(z)

psi1a = @(z) (0.25.*(abs(z).^2) + real(F1(z))); % final interior s'function with Poisson correction (no point vortices)
psi2a = @(z) lambda.*log(abs(z)) + real(F2(z)); % final exterior s'function with far field correction (no point vortices)
psia = @(z) (psi1a(z).*inpolygonc(z,Z) + dsfun(z,zs,lambs) + psi2a(z).*(~inpolygonc(z,Z))); % final combined s'function including pt vortex contribution "dsfun" - see below additional functions

% Velocity u
F1p = @(z) reshape([(0:N).*z(:).^(-1:N-1) -dint./((z(:)-polint).^2)]*c1, size(z)); % analytic function F1'(z)
F2p = @(z) reshape([-(1:N).*(1./z(:)).^(2:N+1) -dext./((z(:)-polext).^2)]*c2, size(z)); % analytic function F2'(z)

u1a = @(z) -1i.*(0.5.*z + conj(F1p(z))); % interior complex velocity with Poisson correction (no point vortices)
u2a = @(z) -1i.*(lambda./conj(z) + conj(F2p(z))); % exterior complex velocity with far field correction (no point vortices)

ua = @(z) u1a(z).*inpolygonc(z,Z) -1i.*nsfun(z,zs,lambs) + u2a(z).*(~inpolygonc(z,Z)); % combined complex velocity including ALL pt vortices "nsfun" - see below additional functions
ua2 = @(z) u1a(z).*inpolygonc(z,Z) -1i.*nsfun2(z,zs,lambs) + u2a(z).*(~inpolygonc(z,Z)); % combined complex velocity including pt vortices EXCLUDING the first one "nsfun2" - see below additional functions (used for calculating velocity imposed on pt vortex 1)

ubdry1 = ua(zseg); ubdry2 = circshift(ubdry1,-1); ubdry = 0.5.*(ubdry1+ubdry2); % complex velocity at each segment boundary point. Following WOZ, the average u_{k+0.5} has been taken - see their (4.1), (4.2) and (4.7R)
uvort = ua2(zs(1)); % complex velocity at first satellite point vortex - by symmetry, all pt vortices have same (magnitude) velocity
end

%% Appendix: Additional functions for satellite point vortices contribution - see Sec. 6 of Overleaf
function dvF = dsfun(z,zs,lambs) % stream function contributions
dvF =  0;
for j=1:size(zs,1)
    dvF = dvF - lambs.*log(abs(z-zs(j)));
end
end

function nvFf = nsfun(z,zs,lambs) % complex velocity contribution - all vorticies
nvFf = 0;
for j=1:size(zs,1)
    nvFf = nvFf - lambs./conj(z-zs(j));
end
end

function nvFf2 = nsfun2(z,zs,lambs) % complex velocity contribution - excluding first vortex
nvFf2 = 0; 
for j=2:size(zs,1)
    nvFf2 = nvFf2 - lambs./conj(z-zs(j));
end
end