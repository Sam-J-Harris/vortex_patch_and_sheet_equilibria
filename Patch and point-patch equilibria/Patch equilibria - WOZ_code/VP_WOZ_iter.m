function [Z, Zseg, OM, AR, PE, RERR] = VP_WOZ_iter(z,zseg,nsym,lambda,omega,itsmax)
% Z = boundary data at each iteration, Zs = segment data at each iteration 
% OM = Omega (angular velocity), AR = area, PE = perimeter area
% RERR = total relative error between iteration shapes
% 
% See appendix A from Harris & McDonald (2025) for details of the WOZ algorithm
%
% Code:
Ra = abs(zseg(1)); Rb = abs(zseg(end)); 
Rab = Ra.^2-Rb^2; M = size(zseg,1)-1; % points span from k=1 to k=M+1
toly = 5*10^(-6); mustar = 0.6; % tolerance in R; relaxation parameter
Rerr = 1; its=1; % terminates at kmax iterations if toly not reached before
Z = z; Zseg = zseg; RERR = []; OM = []; % boundary data, relative error and omega lists
PE = perimeter(polyshape(real(z),imag(z))); AR = polyarea(real(z),imag(z)); % perimeter and area lists

while Rerr>toly && its<=itsmax % terminates once rel error is below toly or once max iterations have been reached
    xcr = real(zseg); xnt = circshift(xcr,-1); xdf = xnt - xcr; % x-difference between points in the segment
    ycr = imag(zseg); ynt = circshift(ycr,-1); ydf = ynt - ycr; % y-difference between points in the segment
    thcr = angle(zseg); thnt = circshift(thcr,-1); thtot = angle(z); % theta lists (theta current, theta next, total theta)
    Rcr = abs(zseg); Rnt = circshift(Rcr,-1); % r lists

    % Step 1: Calculate boundary velocity ubdry - specifically u_(k+0.5)
    [~, ubdry] = VP_WOZ_aaa(z,zseg,lambda); % AAA-LS computation for boundary velocity

    % Step 2: Calculate angular velocity Omega
    Omega = 0; ucr = real(ubdry); vcr = imag(ubdry);
    for k=1:M
        Omega = Omega+(ucr(k).*ydf(k)-vcr(k).*xdf(k));
    end
    Omega = ((2*omega)/Rab).*Omega; OM = [OM Omega]; % find omega, add to big list

    % Step 3: Calculate F^n_(k+0.5) and F^n_(k-1/2)
    Fpl = (ucr.*sin(thnt)-vcr.*cos(thnt)+(Omega/2).*Rnt)./(ucr.*sin(thcr)-vcr.*cos(thcr)+(Omega/2).*Rcr);

    % Step 4: Calculate R^(n+1) = Rcr
    rpts = M-1; H = 0.*ones(rpts,1); % for 2<=k<=M
    H(1) = 0.5.*(Fpl(1)).^(-1).*Ra;
    H(end) = 0.5.*Fpl(M).*Rb;
    A = 0.*ones(rpts,rpts); 
    for k=1:M-1 % from k=2 to k=M
        if k~=1
            A(k,k-1)=-0.5.*(Fpl(k)).^(-1);
        end
            A(k,k)=1;
        if k~=M-1
            A(k,k+1)=-0.5.*(Fpl(k+1));
        end
    end

    Rbr = A\H; Rbr = [Ra; Rbr; Rb];
    Rit = mustar.*Rbr + (1-mustar).*Rcr; % relaxation procedure give R^(n+1)
    
    % Step 5: Update boundary z and boundary segment zs
    Rseg = [Rit; flip(Rit(2:end-1))]; % gives R values on one "fold"
    Rtot = Rseg; 
    for k=1:nsym-1
        Rtot = [Rtot; Rseg]; % add segments
    end
    Rtot = [Rtot; Rtot(1)];
    z = Rtot.*(cos(thtot)+1i*sin(thtot)); Z = [Z z]; % update z, add to big list
    zseg = Rit.*(cos(thcr)+1i*sin(thcr)); Zseg = [Zseg zseg]; % update zseg, add to big list
    
    P = perimeter(polyshape(real(z),imag(z))); PE = [PE P]; % find perimeter, add to big list
    A = polyarea(real(z),imag(z)); lambda = A./(2*pi); AR = [AR A]; % find area, add to big list

    % Step 6: Update Rerr and its
    Rerr = sum(abs(Rit-Rcr)); RERR = [RERR Rerr]; % find new relative error between iterations, add to big list

    if its~=1, fprintf(repmat('\b',1,lineLength)), end
    lineLength= fprintf('WOZ iteration: Step %d of %d completed. Error = %f.\n',its,itsmax,Rerr); % feedback to user on time taken
    its = its+1; % update iterations

%     % Step 7: Plot (OPTIONAL)
%     figure(1+its)
%     plot(real(z),imag(z),'x'); hold on, plot(real(zs),imag(zs),'-x'), hold off
%     daspect([1 1 1]);
end

end