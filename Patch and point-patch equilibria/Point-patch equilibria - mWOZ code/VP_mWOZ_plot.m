function VP_mWOZ_plot(ZAAA,ZSeg,zs,RERR,lambda,lambs,stps,nsym,b,plswt,fignum)
% Plotting function 
% First plot - combined plot of all* iterations (plswt ~= 0)
% Second plot - grid plot of each iteration individually (pltswt == 2)
% Third plot - stream line plot (always plotted)
%
% * up to a tolerance level
%
% Code:
ssz = get(0,'ScreenSize'); mxx = 1.25*abs(zs(1)); mxy = mxx; % capture screen size; set max x and y axis limits
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); % determines if point is in polygon (BOUNDARY IS INCLUDED IN POLYGON) 

if nsym==1, ZAAA = -imag(ZAAA)+1i*real(ZAAA); zs = -imag(zs)+1i*real(zs); end % flips axes for one point - one patch problems
tolyp = 0.5; % tolerance at which it stops plotting
ZApl = ZAAA(:,RERR<tolyp); ZApl = [ZApl, ZAAA(:,end)]; % only plots those with significant difference, plus the last one

% Figure 1 - combined iterations plot
if plswt ~= 0
figure(fignum) % all iterations in one plot
plot(real(ZAAA),-imag(ZAAA),'LineWidth',1.2), hold on, plot(real(zs),-imag(zs),'x','Color','b','MarkerSize',18), hold off
axis([-mxx mxx -mxy mxy]), daspect([1 1 1]);
end

% Figure 2 - all iterations individually in a grid
if plswt == 2
ittot = size(ZAAA,2); itstps = ceil(ittot/stps); m=1; % stepping to get a nice row to column ratio for the grid
colm = ceil(sqrt(itstps)); rows = ceil(itstps/colm); % number of rows and columns
figure(fignum+1)
for j=1:stps:ittot
    subplot(rows,colm,m)
    z = ZAAA(:,j); % plots iteration j in the subplot
    plot(real(z),-imag(z)), hold on, plot(real(zs),-imag(zs),'x'), hold off 
    daspect([1 1 1]), axis([-mxx mxx -mxy mxy]),
    title("Iteration " +num2str(j-1)), m=m+1; % labels the iteration
end
end

% Figure 3 - stream line plot
axlim=3; % meshgrid axis limits (USER INPUT)
levels = linspace(-0.09,0.8,15); % (exterior) contour plot levels (USER INPUT); (-0.1,0.14,10) = 3fold (same sgn); (-0.35,0.04,15) = 3fold (opp sgn); (-0.09,0.8,15) = 4fold (same sgn);
%levels = 50;

z = ZAAA(:,end); zseg = ZSeg(:,end); % takes the final equilibrium shape
[psia,~,uvort] = VP_mWOZ_aaa(z,zseg,zs,lambda,lambs); % AAA-LS computation for s'function and pt vortex velocity - see function for more details
Omega = sign(real(uvort))*abs(uvort)./b; %  condition that pt vortices are fixed on rotating frame
psipl = @(w) psia(w) + (((Omega/2).*abs(w).^2).*(~inpolygonc(w,z)) + ((Omega/2).*abs(w).^2).*(inpolygonc(w,z))); % stream function in the rotating frame

figure(fignum+plswt) % accounts for if other plots have been made for plswt~=0
plot(z,'r','linewidth',4,'Color','r'), hold on, % plots vortex patch boundary  
    
% Interior and exterior s'lines
xs = linspace(-axlim,axlim,500); [xx,yy] = meshgrid(xs,xs); zz = xx+1i*yy;
sfni = psipl(zz); sfne = sfni; sfni(~inpolygonc(zz,z)) = NaN; sfne(inpolygonc(zz,z)) = NaN; 
contour(xs,xs,sfni,5,'linewidth',1,'Color','k'), % INTERIOR stream function contours
contour(xs,xs,sfne,levels,'linewidth',1,'Color','k'), % EXTERIOR stream function contours

axis([-axlim axlim -axlim axlim]), daspect([1 1 1]), % sets the axes 
set(gca,'xdir','reverse','ydir','reverse'), % flip axes vertically
set(gca,'XColor', 'none','YColor','none'), % turn axes off
hold off
end

