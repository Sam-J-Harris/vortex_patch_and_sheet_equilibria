%%  Vortex Patch - modified WOZ method
% Two-domain AAA-LS code for the "modified WOZ" method of Xue et al. 2017
% for steady state solns of N+1 point vortex - vortex patch configs
clc; close all; clear; set(0,'DefaultFigureVisible','on'); warning('off','all'); tic; % clear all, turn warnings off, start timer

%% Parameter Setup (User Input)
Ra = 0.3; % initial circle radius (must be less than b)
nsym = 4; % n-fold symmetry (set = 1 for one pt vortex and patch)
csgn = 2*pi; % circulation of point vortices

M=1000; % M+1= number of points in a segment; 
mustar = 0.6; % relaxation parameter
toly = 5*10^(-6); % tolerance of the iteration (if Rerror < toly, iteration ends)
itsmax = 100; % max iterations (unless steady state reaches toly before itsmax)

%% Create shape data z, segment data zseg and pt vortex data zs
b = 1; % distance of pt vortices from origin (FIX at 1); 
omega = 1; % vorticity of patch (FIX at either 1 or -1);
[z, zseg, zs, lambda, lambs] = VP_mWOZ_shape(Ra,nsym,csgn,M,b,omega); % see shape function for more details

%% AAA-LS + mWOZ iterations
[ZAAA, ZSeg, RERR] = VP_mWOZ_iter(z,zseg,zs,b,nsym,lambda,lambs,M,Ra,itsmax,toly,mustar); % modified WOZ iteration using AAA-LS, gives boundary data list "ZAAA", segment data "ZSeg" and list of relative error between iterations "RERR"
iters = size(ZAAA,2)-1; rtime = toc; % total number of iterations =< itsmax; stop timer
fprintf("Vortex patch mWOZ code (nsym = " + num2str(nsym)+", Ra = " + num2str(Ra)+"): iters = " + num2str(iters) + ", bdry err = " + num2str(RERR(end)) + ", rtime = " +num2str(rtime)+" secs.\n") % output for the user

%% Plot
stps = 1; plswt = 0; % stps=1 to plot all iters; plswt = 0 (only show s'line plots), = 1 (show combined iters and s'line plot), = 2 (show combined + individual iters and s'line plot)
VP_mWOZ_plot(ZAAA,ZSeg,zs,RERR,lambda,lambs,stps,nsym,b,plswt,1) % plotting function - see function for more details