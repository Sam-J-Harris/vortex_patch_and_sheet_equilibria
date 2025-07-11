function [REr, REo, REa, REp] = VP_WOZ_err(nsym,shswt,OM,AR,PE,RERR)
% = finds the relative errors of: the iterations, angular velocity (omega), area and perimeter
Otab = [[0.33310 0.32843 0.31996 0.30948 0.30122]; [0.37460 0.37153 0.36620 0.35938 0.35395]; [0.39945 0.39705 0.39313 0.38815 0.38425]; [0.41596 0.41395 0.41087 0.4070 0.40407]]; % omega values from WOZ Table III
Atab = [[3.2967 3.7967 4.2382 4.5909 4.7291]; [3.2947 3.5684 3.8046 3.9870 4.0546]; [3.2927 3.4726 3.6263 3.7426 3.7842]; [3.2909 3.4193 3.5284 3.6097 3.6375]]; % area values from WOZ Table III
Ptab = [[6.4440 7.0463 7.7082 8.4021 9.0652]; [6.4488 6.8154 7.2310 7.6811 8.1380]; [6.4555 6.7259 7.0389 7.3859 7.7534]; [6.4641 6.6827 6.9388 7.2271 7.5415]]; % perimeter values from WOZ Table III

Ome = Otab(max(nsym-2,1),shswt); Ae = Atab(max(nsym-2,1),shswt); Pe = Ptab(max(nsym-2,1),shswt); % exact values from WOZ
Omf = -OM(end); Af = AR(end); Pf = PE(end); % values from AAA-LS--WOZ code (omega is negative)
REr = RERR(end); REo = reftn(Ome,Omf); REa = reftn(Ae,Af); REp = reftn(Pe,Pf); % relative errors for plots, omega value, area, perimeter
end

%% Appendix: Relative Error function
function EC = reftn(uex, uap)
EC = abs(uex-uap)./abs(uex); % relerr = abs(exact - approx)/abs(exact). Difficulty if either exact=0 or approx = 0.
end