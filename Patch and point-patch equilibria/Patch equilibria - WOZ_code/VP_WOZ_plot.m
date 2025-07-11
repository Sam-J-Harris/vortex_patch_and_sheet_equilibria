function VP_WOZ_plot(ZAAA,ZAseg,RERR,omega,lambda,stps,plswt,fignum)
% displays plot of vortex patch and interior/exterior stream lines
ssz = get(0,'ScreenSize'); set(gca,'FontSize',20);
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); % determines if point is in polygon (BOUNDARY IS INCLUDED IN POLYGON) 

tolyp = 0.5; % tolerance at which it stops plotting
ZApl = ZAAA(:,RERR<tolyp); ZApl = [ZApl, ZAAA(:,end)]; % only plots patch boundaries with significant differences, plus the last one

% Figure 1 - combined iters and first and last iter (optional)
if plswt ~= 0
figure(fignum)
subplot(3,7,[1,2,3,8,9,10,15,16,17]);
plot(real(ZApl),imag(ZApl),'LineWidth',1.2), 
scl = 1.25; mxx = scl.*max(max(abs(real(ZApl)))); mxy = scl.*max(max(abs(imag(ZApl)))); 
daspect([1 1 1]), axis([-mxx mxx -mxy mxy]), set(gca,'FontSize',20),
set(gca,'XColor', 'none','YColor','none'), 

for j=1:2
    if j==1
        pos = [11,12,18,19]; z = ZAAA(:,1); word = "1";
    else     
        pos = [13,14,20,21]; z = ZAAA(:,end); word = num2str(size(ZAAA,2)-1);
    end
    subplot(3,7,pos)
    x = real(z); y = imag(z); mxx = 1.1*max(abs(x)); mxy = 1.1*max(abs(y));
    plot(x,y,'LineWidth',2), daspect([1 1 1]), axis([-mxx mxx -mxy mxy]), set(gca,'FontSize',20), set(gca,'XColor', 'none','YColor','none'), 
    title("Step " +word,'FontSize',20),
end
hold off
end

% Figure 2 - all iterations (optional)
if plswt == 2
ittot = size(ZAAA,2); itstps = ceil(ittot/stps); m=1;
colm = ceil(sqrt(itstps)); rows = ceil(itstps/colm); 
figure(fignum+1)
for j=1:stps:ittot
    subplot(rows,colm,m)
    z = ZAAA(:,j); x = real(z); y = imag(z); mxx = 1.1*max(abs(x)); mxy = 1.1*max(abs(y));
    plot(real(z),imag(z)), daspect([1 1 1]), axis([-mxx mxx -mxy mxy]),
    title("Iteration " +num2str(j-1)), m=m+1;
end
set(gcf, 'Position',  [ssz(3)/2, ssz(4)/4.5, ssz(3)/2.5, ssz(3)/3])
end

% Figure 3 - stream line plot of final iteration
z = ZAAA(:,end); zseg = ZAseg(:,end);
Ra = abs(zseg(1)); Rb = abs(zseg(end));
Rab = Ra.^2-Rb^2; M = size(zseg,1)-1; % points span from k=1 to k=M+1

xcr = real(zseg); xnt = circshift(xcr,-1); xdf = xnt - xcr;
ycr = imag(zseg); ynt = circshift(ycr,-1); ydf = ynt - ycr;

[psia, ubdry] = VP_WOZ_aaa(z,zseg,lambda); % Step 1: Calculate boundary velocity ubdry - specifically u_(k+0.5)

Omega = 0; ucr = real(ubdry); vcr = imag(ubdry);
for k=1:M, Omega = Omega+(ucr(k).*ydf(k)-vcr(k).*xdf(k)); end, Omega = ((2*omega)/Rab).*Omega; % Step 2: Calculate angular velocity Omega
psipl = @(w) psia(w) + ((Omega/2).*abs(w).^2).*(~inpolygonc(w,z)) + ((Omega/2).*abs(w).^2).*(inpolygonc(w,z)); % stream function used for plotting

figure(fignum+plswt)
axlim=2.5;  
plot(z,'r','linewidth',4,'Color','r'), hold on, % plots vortex patch boundary

% Interior and exterior s'lines
xs = linspace(-axlim,axlim,500); [xx,yy] = meshgrid(xs,xs); zz = xx+1i*yy;
sfni = psipl(zz); sfne = sfni; sfni(~inpolygonc(zz,z)) = NaN; sfne(inpolygonc(zz,z)) = NaN; 
levels = linspace(-0.3,-0.05,15); % -0.3,-0.05,15: for 3-fold; -0.4,-0.125,15: for 4-fold
contour(xs,xs,sfni,5,'linewidth',1,'Color','k'), % INTERIOR stream function contours
contour(xs,xs,sfne, levels,'linewidth',1,'Color','k'), % EXTERIOR stream function contours

axis([-axlim axlim -axlim axlim]), axis square, set(gca,'FontSize',20),
set(gca,'XColor', 'none','YColor','none'), %set(gca, 'color', 'none');
hold off
end