clear; clc; close all;

%This code generates the initial conditions and grid for the CFD solver. 
%This part was transferred to Matlab because it's more readable. 
%Also an external interface for this section is useful because it is not
%hard-coded!

dt=1e-4;
dx=5/150;
dy=5/150;
nx=150;
ny=150;
Re=200;

last_t=1e6;
t_decimation=10;
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
uRightKind=Dirichlet;
uTopKind=Dirichlet;
uBottomKind=Dirichlet;

vLeftKind=Dirichlet;
vRightKind=Dirichlet;
vTopKind=Dirichlet;
vBottomKind=Outflow;

%Domain
u=ones(ny,nx-1)*0;
v=ones(ny-1,nx)*0;
p=ones(ny,nx)*0;

%Saves
save('Initialization.mat');
disp('saved!')