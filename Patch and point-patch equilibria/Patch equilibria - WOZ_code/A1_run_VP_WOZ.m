%%  Vortex Patch - WOZ method
% Finds n-fold symmetric vortex patch using the WOZ numerical algorithm.
% Requires the chebfun package: download from https://www.chebfun.org/
%
% Code:
clc; close all; clear; set(0,'DefaultFigureVisible','on'); warning('off','all'); tic; % start timer

%% Parameter Setup (User Input)
nsym = 3; % n-fold symmetry - either 3, 4, 5 or 6
shswt = 3; % choice of Ra - see shape function
M=120; % number of pts in a segment
itsmax = 500;  % max iterations
Pts = M*(2*nsym)+1; % total number of boundary pts; 
[z, zs, Ra, A] = VP_WOZ_shape(nsym,shswt,M,Pts); % boundary data, segment data and area of patch
omega = 1; % vorticity (in rotations/second [1/s]);
lambda = A/(2*pi); % (initial) dimless constant

%% AAA-LS Iterations
[ZAAA, ZAseg, OM, AR, PE, RERR] = VP_WOZ_iter(z,zs,nsym,lambda,omega,itsmax); % WOZ iteration using AAA-LS, gives boundary data, Omega (angular velocity), area, perimeter, and relative error between iterations
rtime = toc; fprintf("Vortex patch WOZ code for nsym = " + num2str(nsym)+" and Ra = " + num2str(Ra)+": steps = " + num2str(size(ZAAA,2)-1)+", runtime = " +num2str(rtime)+" seconds.\n")

%% Comparison with WOZ results
[REr, REo, REa, REp] = VP_WOZ_err(nsym,shswt,OM,AR,PE,RERR); % finds the relative errors of: the iterations, angular velocity (omega), area and perimeter
fprintf("Relative errors: iteration = "+num2str(REr)+", omega = "+num2str(REo)+", area = "+num2str(REa)+", perimeter = "+num2str(REp)+".\n")

%% Plot
stps = 50; plswt = 1; % number of iterations to be plotted; plswt = 0 (only show s'line plots), = 1 (show combined iters and s'line plot), = 2 (show combined + individual iters and s'line plot)
VP_WOZ_plot(ZAAA,ZAseg,RERR,omega,lambda,stps,plswt,1); % plotting function

