clear; clc; close all;

%This code generates the initial conditions and grid for the CFD solver. 
%This part was transferred to Matlab because it's more readable. 
%Also an external interface for this section is useful because it is not
%hard-coded!

nRef=1;%Ratio of Re/100

dt=5e-3/nRef;
nx=300*nRef;
ny=260*nRef;
dx=15/nx;
dy=13/ny;
Re=100*nRef;

uTest=3;
CellRe=Re*dx*uTest %u is considered sqrt(2)
CFL=dt*uTest/dx

last_t=50e3;
t_decimation=100;
reltol=1e-6;

%Boundary conditions
uLeft=ones(1,ny)*0;
uRight=ones(1,ny)*0;
uTop=ones(1,nx-1)*0;
uBottom=ones(1,nx-1)*0;

vLeft=ones(1,ny-1)*0;
vRight=ones(1,ny-1)*0;
vTop=ones(1,nx)*0;
vBottom=ones(1,nx)*0;

%Boundary kinds
Dirichlet=0;
Outflow=1;

uLeftKind=Dirichlet;
uRightKind=Outflow;
uTopKind=Dirichlet;
uBottomKind=Dirichlet;

vLeftKind=Dirichlet;
vRightKind=Dirichlet;
vTopKind=Dirichlet;
vBottomKind=Dirichlet;

%Domain
% [x, y]=ndgrid(((1:nx)-1)/nx,((1:ny)-1)/ny);
% d=sqrt((x-0.5).^2+(y-0.5).^2);
% e=exp(-(d.^2)*10);
% theta=atan2(y-0.5, x-0.5);

u=randn(ny,nx-1)*0.05; %exp(-(d.^2)*10).*cos(theta);
v=randn(ny-1,nx)*0;%-exp(-(d.^2)*10).*sin(theta);
p=ones(ny,nx)*0;

%Saves
save('Initialization.mat');
disp('saved!')