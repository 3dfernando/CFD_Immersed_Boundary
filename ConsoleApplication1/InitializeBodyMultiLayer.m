clear; clc; close all;

%%
%This code creates a file 'Body.mat' that contains arrays with the
%coordinates of the points of the bodies in the flow.
%The input is a DXF file that contains the lines that constitute the body
%boundaries.
%For multiple bodies, use multiple layers in the DXF.

%Velocity and trajectories need to be coded here in the next section of the
%code!

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
        Layers{ns}=CirclesLayers{i};
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
        Layers{ns}=ArcsLayers{i};
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
        Layers{ns}=LinesLayers{i};
    end 
end


for k=1:size(c_Poly,1)
    %Generates Polylines
    Poly=c_Poly{k,1};
    for i=1:(size(Poly,1)-1)
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
            Layers{ns}=c_Poly{k,2};
        end 
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
        Layers{i}=[];
    end
end

%Removes points too close to the bounds
Lx=dx*nx;
Ly=dy*ny;
padX=2*dx;
padY=2*dy;

yP(xP<padX)=nan; Layers(xP<padX)=cell(1,sum(xP<padX)); xP(xP<padX)=nan; %Left bound
yP(xP>(Lx-padX))=nan; Layers(xP>(Lx-padX))=cell(1,sum(xP>(Lx-padX))); xP(xP>(Lx-padX))=nan; %Right bound
xP(yP<padY)=nan; Layers(yP<padY)=cell(1,sum(yP<padY)); yP(yP<padY)=nan; %Bottom bound
xP(yP>(Ly-padY))=nan; Layers(yP>(Ly-padY))=cell(1,sum(yP>(Ly-padY))); yP(yP>(Ly-padY))=nan; %Bottom bound

%Actually removes the points here (nans)
yP=yP(~isnan(xP));
xP=xP(~isnan(xP));
Layers = Layers(~cellfun('isempty',Layers));
ns=length(xP);

UniqueLayers=unique(Layers); %For object selection



%%
%=================TRAJECTORY CODING=====================
%Creates a trajectory for the points
time=dt*((1:last_t)-1);
time=repmat(time,ns,1);

xTraj=time*0;%Inits
yTraj=time*0;%Inits

%Each layer can have its own trajectory
%=======Layer 0========
i=1; %First object
currentLayer = strcmp(Layers,UniqueLayers{i});
xTraj(currentLayer,:)=0;
yTraj(currentLayer,:)=0;

%=======Layer 1========
i=2; %Second object
currentLayer = strcmp(Layers,UniqueLayers{i});
xTraj(currentLayer,:)=0.5*sin(2*pi*0.05*time(currentLayer,:)); %Moves piston back and forth
yTraj(currentLayer,:)=0;


%saves the data as a mat file
save('Body.mat','xP','yP','xC','yC','xTraj','yTraj', 'ns', 'ds');
disp('saved!')



%%
%Extras!
%Calculates vx and vy maxima (just for my information)
vx=(xTraj(:,2:end)-xTraj(:,1:end-1))/dt;
vy=(yTraj(:,2:end)-yTraj(:,1:end-1))/dt;

vxMax=max(vx(:))
vyMax=max(vy(:))

%=======================================================================
%Plots the points
colors={'k.','b.','r.','g.','k*','b*','r*','g*'};

figure; 
for t_idx=1:100:last_t
    cla; hold on;
    for i=1:length(UniqueLayers)
        currentLayer = strcmp(Layers,UniqueLayers{i});
        plot(xP(currentLayer) + xTraj(currentLayer,t_idx),yP(currentLayer) + yTraj(currentLayer,t_idx),colors{i});
        daspect([1 1 1]);
    end
    plot([0 Lx Lx 0 0], [0 0 Ly Ly 0],'b:');
    plot(xC,yC,'b*');
    hold off;
    drawnow;
    
    pause(0.1)
end



