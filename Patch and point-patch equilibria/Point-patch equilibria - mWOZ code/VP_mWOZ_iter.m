function [Z, ZSeg, RERR] = VP_mWOZ_iter(z,zseg,zs,b,nsym,lambda,lambs,M,Ra,itsmax,toly,mustar)
% Performs modified WOZ iteration of Xue et al. (2017)
% Iterates towards steady state of point vortex - vortex patch config
%
% See Sec. 6 and Appendix A from Harris & McDonald (2025) for details of the mWOZ algorithm
%
% Code:
Rerr = 1; its=1; % initial error in R; iteration counter (if its>itsmax, iteration ends)
Z = z; ZSeg = zseg; RERR = []; % boundary data, relative error list

while Rerr>toly && its<=itsmax
    thcr = angle(zseg); thnt = circshift(thcr,-1); thtot = angle(z); % list of next theta values in the segment; theta values of entire contour
    Rcr = abs(zseg); Rnt = circshift(Rcr,-1); % list of next R values in the segment

    %%% STEPS FROM WOZ SECTION 4.4 %%%
    % Step 1: Calculate boundary velocity ubdry - specifically u_(k+0.5)
    [~,ubdry, uvort] = VP_mWOZ_aaa(z,zseg,zs,lambda,lambs); % AAA-LS computation for bdry and pt vortex velocity - see function for more details
    ucr = real(ubdry); vcr = imag(ubdry); % extract velocities u and v from complex velocity ubdry = u + iv

    % Step 2: Calculate angular velocity Omega (main difference from standard WOZ)
    Omega = sign(real(uvort))*abs(uvort)./b; %  condition that pt vortices are fixed on rotating frame - uvort = complex velocity of point vortices

    % Step 3: Calculate F^n_(k+0.5) and F^n_(k-1/2)
    Fpl = (ucr.*sin(thnt)-vcr.*cos(thnt)+(Omega/2).*Rnt)./(ucr.*sin(thcr)-vcr.*cos(thcr)+(Omega/2).*Rcr); % equation (4.7R) from WOZ

    % Step 4: Calculate R^(n+1) = Rcr (new current R) - equation (4.8) from WOZ
    rpts = M-1; H = 0.*ones(rpts,1); % for 2<=k<=M
    H(1) = 0.5.*(Fpl(1)).^(-1).*Ra;
    H(end) = 0.5.*Fpl(M).*Ra; % instead of (4.8) = 0, put terms involving R_1 = Ra and R_M+1 = Rb(=Ra for mWOZ) on the right hand side, so (4.8) = H
    A = 0.*ones(rpts,rpts); % AR = H where A = matrix of coeffs (involving Fpl), R = vector of R's (to be found), H = right hand side, as line above
    for k=1:M-1 % from k=2 to k=M
        if k~=1
            A(k,k-1)=-0.5.*(Fpl(k)).^(-1);
        end
            A(k,k)=1;
        if k~=M-1
            A(k,k+1)=-0.5.*(Fpl(k+1));
        end
    end % correct coefficients of A for corresponding R_k

    Rbr = A\H; Rbr = [Ra; Rbr; Ra]; % calculates R bar as in (4.8), noting R_1 = R_M+1 = Ra
    Rit = mustar.*Rbr + (1-mustar).*Rcr; % relaxation procedure give R^(n+1)
    
    % Step 5: Update boundary z and boundary segment zs
    Rtot = Rit(1:end-1); % gives R values on one "fold" (remove M+1 point: only M points on the segment)
    for k=1:nsym-1
        Rtot = [Rtot; Rit(1:end-1)]; % add segments together to get total R values
    end
    Rtot = [Rtot; Rtot(1)]; % close the curve
    z = Rtot.*(cos(thtot)+1i*sin(thtot)); zseg = Rit.*(cos(thcr)+1i*sin(thcr)); % update z (total curve coordinates) and zseg (segment coordinates)
    Z = [Z z]; ZSeg = [ZSeg zseg]; % update boundary data and segment data lists
    A = polyarea(real(z),imag(z)); lambda = A./(2*pi*b^2); % update lambda value

    % Step 6: Update Rerr and its
    Rerr = sum(abs(Rit-Rcr)); RERR = [RERR Rerr]; % error function - equation (4.17) of WOZ. Add this error to the list of errors

    if its~=1, fprintf(repmat('\b',1,lineLength)), end
    lineLength= fprintf('WOZ iteration: Step %d of %d completed. Error = %f.\n',its,itsmax,Rerr); %feedback to user on time taken
    its = its+1; % update iterations

    % % Step 7: Plot (OPTIONAL)
    % figure(1+its)
    % plot(real(z),imag(z)); hold on, plot(real(zs),imag(zs),'x'), hold off
    % daspect([1 1 1]);
end

end