function [psia, ubdry] = VP_WOZ_aaa(Z,Zs,lambda)
% psia = stream function, ubdry = velocity on the patch boundary
%
% Code:
% AAA-LS numerical method
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); % determines if point is in polygon (BOUNDARY IS INCLUDED IN POLYGON)

% Normal vector
dZ = circshift(Z(1:end-1),1) - circshift(Z(1:end-1),-1); v = 1i.*(dZ)./abs(dZ); 
nvec = v./abs(v); nvec(end+1) = nvec(1); % unit normal vector

% Far field condition into Dirichlet and Neumann BCs
hd = @(z) lambda.*log(abs(z)) - 0.25.*(abs(z).^2); % Dirichlet BC 
hn = @(z) lambda.*real(nvec./z) - 0.5.*(real(nvec.*conj(z))); % Neumann BC
Hd = hd(Z); Hn = hn(Z); bigH = [Hd; Hn]; 

% Global AAA poles %%%%%%%% 
rtoly = 1e-8; ztoly = 1e-2;
[~,poltot1,restot1] = aaa(Hd,Z,'cleanup',1,'toly',rtoly); % poles for Dirichlet BC
[~,poltot2,restot2] = aaa(Hn,Z,'cleanup',1,'toly',rtoly); % poles for Neumann BC

poltot = [poltot1; poltot2].'; restot = [restot1; restot2].'; % all poles and resolutions

% Pole control - manual elimination of Froissart doublets  %%%%%%%%
poltot = poltot(abs(restot)>rtoly); % poles of small residue < rtoly eliminated
D = min(abs(Z-poltot),[],1); D(D<ztoly*max(abs(Z))) = 0; D(D~=0) =1; % finds if min dist. b/ween poles and bdry is below toly
poltot=poltot.*D; poltot(abs(poltot)==0)=[]; % poles close to the boundary < ztoly eliminated

polint = (poltot(~inpolygonc(poltot,Z))); % exterior poles (used in interior pb)
polext = (poltot(inpolygonc(poltot,Z))); % interior poles (used in exterior pb)

% % Plot of poles %%%%%%%%
% figure(3) % uncomment if you want to see figure of the poles
% zplot = Z(1:end); xp = real(zplot); yp = imag(zplot);
% plot(xp,yp,'k'), hold on, plot(polint,'.r'); plot(polext,'.b'), hold off, 
% axis square, daspect([1 1 1]); axis([-10 10 -10 10]);

% Interior and exterior poles with min distance to bdry
if size(polint,1)==0 && size(polint,2)==0 % interior poles
    polint=0; dint=0;
else
    dint = min(abs((Z)-polint),[],1); % d=distance between pole and closest point on polygon
end

if size(polext,1)==0 && size(polext,2)==0 % exterior poles
    polext=0; dext=0;
else
    dext = min(abs((Z)-polext),[],1); 
end

% Least-Squares %%%%%%%% 
N = 20; N2int = size(polint,2);
P1d = Z.^(0:N); % interior Dirichlet Runge bases
P2d = (1./Z).^(1:N); % exterior Dirichlet Runge bases

P1n = (0:N).*Z.^(-1:N-1).*nvec; % interior Neumann Runge bases
P2n = -(1:N).*Z.^(-(2:N+1)).*nvec; % exterior Neumann Runge bases

Q1d = dint./((Z-polint)); % interior Dirichlet Newman bases
Q2d = dext./((Z-polext)); % exterior Dirichlet Newman bases

Q1n = -dint./((Z-polint).^2).*nvec; % interior Neumann Newman bases
Q2n = -dext./((Z-polext).^2).*nvec; % exterior Neumann Newman bases
      
A1d = [real(P1d) real(Q1d) -imag(P1d) -imag(Q1d)];
A2d = [real(P2d) real(Q2d) -imag(P2d) -imag(Q2d)]; 

A1n = [real(P1n) real(Q1n) -imag(P1n) -imag(Q1n)];
A2n = [real(P2n) real(Q2n) -imag(P2n) -imag(Q2n)];

bigA = [A1d A2d; A1n A2n]; testc = bigA\bigH; % least-squares operation

c1 = testc(1:(2*(N+1)+2*N2int),:); c2 = -testc((2*(N+1)+2*N2int)+1:end,:);
c1 = reshape(c1,[],2)*[1;1i]; c2 = reshape(c2,[],2)*[1;1i]; % complex coefficients

% Stream function psi
F1 = @(z) reshape([z(:).^(0:N) dint./(z(:)-polint)]*c1, size(z)); %analytic function F1(z)
F2 = @(z) reshape([(1./z(:)).^(1:N) dext./(z(:)-polext)]*c2, size(z)); %analytic function F2(z)

psi1a = @(z) (0.25.*(abs(z).^2) + real(F1(z))); % final stream function includes additional terms (Poisson and far-field corrections)
psi2a = @(z) lambda.*log(abs(z)) + real(F2(z));
psia = @(z) (psi1a(z).*inpolygonc(z,Z) + psi2a(z).*(~inpolygonc(z,Z))); % stream function 

% Velocity u
F1p = @(z) reshape([(0:N).*z(:).^(-1:N-1) -dint./((z(:)-polint).^2)]*c1, size(z)); % analytic function F1'(z)
F2p = @(z) reshape([-(1:N).*(1./z(:)).^(2:N+1) -dext./((z(:)-polext).^2)]*c2, size(z)); % analytic function F2'(z)

u1a = @(z) -1i.*(0.5.*z + conj(F1p(z))); % interior potential V1
u2a = @(z) -1i.*(lambda./conj(z) + conj(F2p(z))); % exterior potential V2
ua = @(z) u1a(z).*inpolygonc(z,Z) + u2a(z).*(~inpolygonc(z,Z)); % velocity

ubdry1 = ua(Zs); ubdry2 = circshift(ubdry1,-1); ubdry = 0.5.*(ubdry1+ubdry2); % circshift is -1 to bring next data point to correct place
end