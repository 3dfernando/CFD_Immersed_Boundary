clear; clc; close all;

load('Initialization.mat');

[f,p]=uigetfile('*.dxf','Please select a drawing file');
[c_Line,c_Poly,c_Cir,c_Arc,c_Poi] = readDXF(fullfile(p,f));

%Gets points coordinates
Circles=cell2mat(c_Cir(:,1));
Lines=cell2mat(c_Line(:,1));
Polylines=cell2mat(c_Poly(:,1));
Arcs=cell2mat(c_Arc(:,1));
Points=cell2mat(c_Poi(:,1));

%Gets layers [LAYER 0 IS VALID OBJECTS, OTHER LAYERS ARE TRAJECTORIES
CirclesLayers=c_Cir(:,2);
LinesLayers=c_Line(:,2);
PolylinesLayers=c_Poly(:,2);
ArcsLayers=c_Line(:,2);
PointsLayers=c_Poi(:,2);


%Generates points according to ds
ds=min(dx,dy);
ns=0;

for i=1:size(Circles,1)
    %Generates circles
    r=Circles(i, 3);
    nTheta=ceil(2*pi*r/ds);
    theta=linspace(0, 2*pi, nTheta);
    theta=theta(1:end-1);
    
    xCenter=Circles(i, 1);
    yCenter=Circles(i, 2);
    for j=1:length(theta)
        ns=ns+1;
        xP(ns)=xCenter + r*cos(theta(j));
        yP(ns)=yCenter + r*sin(theta(j));
    end 
end

for i=1:size(Arcs,1)
    %Generates arcs
    r=Arcs(i, 3);
    angleStart=deg2rad(Arcs(i, 4));
    angleEnd=deg2rad(Arcs(i, 5));
    
    if angleStart>angleEnd
        angleEnd = angleEnd + 2*pi;
    end
    
    nTheta=ceil(r*abs(angleEnd-angleStart)/ds);
    theta=linspace(angleStart, angleEnd, nTheta);
    
    xCenter=Arcs(i, 1);
    yCenter=Arcs(i, 2);
    for j=1:length(theta)
        ns=ns+1;
        xP(ns)=xCenter + r*cos(theta(j));
        yP(ns)=yCenter + r*sin(theta(j));
    end 
end


for i=1:size(Lines,1)
    %Generates lines
    x0=Lines(i,1);
    y0=Lines(i,2);
    
    x1=Lines(i,4);
    y1=Lines(i,5);
    
    lineLength=sqrt((x1-x0)^2 + (y1-y0)^2);
    nSegments=ceil(lineLength/ds);
    segments=linspace(0, 1, nSegments);
    
    for j=1:length(segments)
        ns=ns+1;
        xP(ns)=x0 + (x1-x0)*segments(j);
        yP(ns)=y0 + (y1-y0)*segments(j);
    end 
end


for i=1:(size(Polylines,1)-1)
    %Generates Polylines
    x0=Polylines(i,1);
    y0=Polylines(i,2);
    
    x1=Polylines(i+1,1);
    y1=Polylines(i+1,2);
    
    lineLength=sqrt((x1-x0)^2 + (y1-y0)^2);
    nSegments=ceil(lineLength/ds);
    segments=linspace(0, 1, nSegments);
    
    for j=1:length(segments)
        ns=ns+1;
        xP(ns)=x0 + (x1-x0)*segments(j);
        yP(ns)=y0 + (y1-y0)*segments(j);
    end 
end

%First point will be used as model center
if(~isempty(Points))
    xC=Points(1,1);
    yC=Points(1,2);
else
    xC=0;
    yC=0;
end


%Removes points too close together
threshold=ds/2;
D = pdist2([xP.' yP.'],[xP.' yP.'],'euclidean');
tooClose=D<threshold;
tooClose=tril(tooClose,-1);

for i=1:size(tooClose,1)
    %Deletes each of the repeated points
    if sum(tooClose(i,:))~=0
        %Means this is a repeated point
        tooClose(i,:)=0;
        tooClose(:,i)=0;
        xP(i)=nan;
        yP(i)=nan;
    end
end

%Removes points too close to the bounds
Lx=dx*nx;
Ly=dy*ny;
padX=2*dx;
padY=2*dy;

yP(xP<padX)=nan; xP(xP<padX)=nan; %Left bound
yP(xP>(Lx-padX))=nan; xP(xP>(Lx-padX))=nan; %Right bound
xP(yP<padY)=nan; yP(yP<padY)=nan; %Bottom bound
xP(yP>(Ly-padY))=nan; yP(yP>(Ly-padY))=nan; %Bottom bound

%Actually removes the points here (nans)
yP=yP(~isnan(xP));
xP=xP(~isnan(xP));
ns=length(xP);

%Plots the points
plot(xP,yP,'k.');
hold on;
plot([0 Lx Lx 0 0], [0 0 Ly Ly 0],'b:');
plot(xC,yC,'b*');
daspect([1 1 1]);


%Creates a trajectory for the points
time=dt*((1:last_t)-1);
xTraj=time*0;
yTraj=time*0;

f=0.5;
thetaTraj=time*(2*pi*f); %Theta in rad
thetaTraj(time>(0.125/f))=(0.125/f)*(2*pi*f); %Stops in the middle
thetaTraj=-thetaTraj;


%saves the data as a mat file
save('Body.mat','xP','yP','xC','yC','xTraj','yTraj','thetaTraj', 'ns', 'ds');
disp('saved!')


