%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 1 September 2015.

%%

addpath ../../grid/
addpath ../../basis/
addpath ../../particleGen/
addpath ../../constitutiveModels/
addpath ../../util/
addpath ../../geoMesh/
addpath ../../externals/PolyMesher/
addpath ../../externals/
addpath ../../postProcessing/

%%
clc
colordef white


opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

rho    = 2.5;

tic;

%% Computational grid (all length in cm)

l = 4;
w = 4;

noX0      = 42;        % number of elements along X direction
noY0      = 42;        % number of elements along Y direction
ghostCell = 0;

[grid]    = buildGrid2D(l,w,noX0,noY0, ghostCell);

node      = grid.node;
element   = grid.element;
deltax    = grid.deltax;
deltay    = grid.deltay;
elemCount = grid.elemCount;
nodeCount = grid.nodeCount;
numx2     = grid.numx;
numy2     = grid.numy;
Vc        = deltax * deltay;

%% generate material points

ppc           = [3 3];
circle.center = [2 2];
circle.radius = 1;  % cm

[res]          = generateMPForCircle(circle,ppc,grid);

pCount         = size(res.position,1);
body1.volume   = res.volume;
body1.volume0  = res.volume;
body1.mass     = res.volume*rho;
body1.coord    = res.position;
body1.deform   = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress   = zeros(pCount,3);                % stress
body1.strain   = zeros(pCount,3);                % strain
body1.velo     = zeros(pCount,2);                % velocity
body1.pressure = zeros(pCount,1);                % pressure
body1.density  = rho*ones(pCount,1);             % density (change in time)


% put all bodies in one variable
bodies    = cell(1,1);
bodies{1} = body1;
bodyCount = length(bodies);

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

for ib=1:length(bodies)
  body      = bodies{ib};
  elems     = ones(length(body.volume),1);
  
  for ip=1:length(body.volume)
    x = body.coord(ip,1); y = body.coord(ip,2);
    e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
    elems(ip) = e;
  end
  
  bodies{ib}.elements = unique(elems);
  bodies{ib}.nodes    = unique(element(bodies{ib}.elements,:));
  mpoints = cell(elemCount,1);
  for ie=1:elemCount
    id  = find(elems==ie);
    mpoints{ie}=id;
  end
  
  bodies{ib}.mpoints  = mpoints;
end

% boundary nodes

%% node quantities

nmassS    = zeros(nodeCount,1);  % nodal mass vector of the system
nmomentumS= zeros(nodeCount,2);  % nodal momentum vector of the system
niforceS  = zeros(nodeCount,2);  % nodal internal force of the system
nmass     = zeros(nodeCount,1);  % nodal mass vector of each body
nmomentum = zeros(nodeCount,2);  % nodal momentum vector of each body
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector
nvelo     = zeros(nodeCount,2*(bodyCount+1));  % nodal velocities (body1,body2,center of mass)
nvelo0    = nvelo;
nacce     = zeros(nodeCount,2*bodyCount);



%% check grid normals in a multiple body case

figure(100)
hold on
coords=bodies{1}.coord;
plot_mesh(node,element,'Q4','k-',1.);
plot(coords(:,1),coords(:,2),'k.','markersize',10);

tol = 1e-15;
for ib=1:bodyCount
  body      = bodies{ib};
  %nodes     = body.nodes;
  [res]     = computeGridNormalCurvature(grid,body);
  bodies{ib}.normals    = res.gNormals;
  bodies{ib}.curvatures = res.gCurvatures;
  normals    = res.gNormals;
  density    = res.cDensities;
  curvature  = res.cCurvatures;
  gcurvature = res.gCurvatures;
  for i=1:length(node)
    nid = i;
    xI  = node(nid,:);
    nI  = normals(nid,:);
   % if norm(nI) > tol, nI = nI / norm(nI); end
    le = 0.01;    
    quiver(xI(1), xI(2), le*nI(1), le*nI(2),4,'LineWidth',2,...
        'Color','r','MaxHeadSize',10);      
  end
end
axis off
title('Un-normalized grid normals', 'FontSize', 18)
%%
figure
hold on
plot_field(node,element,'Q4',sqrt(sum(normals.^2,2)));
%plot_field(node,element,'Q4',curvature);
colormap(jet(10));
colorbar
h=colorbar;
set(h,'fontsize',14);
axis off
title('Length of grid normals', 'FontSize', 18)

%%
figure
hold on
plot_field(node,element,'Q4',gcurvature);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
colormap jet
colorbar
h=colorbar;
set(h,'fontsize',14);
axis off
title('Grid curvature', 'FontSize', 18)

figure
hold on
[X,Y] = meshgrid(0:l/noX0:l,0:l/noY0:l);
C     = reshape(-gcurvature,noX0+1,noY0+1);
surf(X,Y,C)
colormap jet
colorbar
h=colorbar;
set(h,'fontsize',14);
axis on
view(3)
title('Grid curvature', 'FontSize', 18)

%%

figure
hold on
patch('Faces',element,'Vertices',node,'FaceVertexCData',density,'FaceColor','flat')
plot(coords(:,1),coords(:,2),'k.','markersize',10);
%caxis([0.0 1]); 
colormap(jet(10)); % <- change color reps...
h=colorbar;
set(h,'fontsize',14);
axis equal
axis off
title('Cell-centered density', 'FontSize', 18);

figure
hold on
patch('Faces',element,'Vertices',node,'FaceVertexCData',res.scDensities,'FaceColor','flat')
plot(coords(:,1),coords(:,2),'k.','markersize',10);
%caxis([0.0 1]); 
colormap(jet(10)); % <- change color reps...
h=colorbar;
set(h,'fontsize',14);
axis equal
axis off
title('Smoothed cell-centered density', 'FontSize', 18);

%%

figure
hold on
patch('Faces',element,'Vertices',node,'FaceVertexCData',-curvature,'FaceColor','flat')
plot(coords(:,1),coords(:,2),'k.','markersize',10);
%caxis([0.0 1]); 
colormap(jet(10)); % <- change color reps...
h=colorbar;
set(h,'fontsize',14);
axis equal
axis off
title('Cell-centered curvature', 'FontSize', 18);

%%

figure 
hold on
plot_mesh(node,element,'Q4','k-',1.);
plot(coords(:,1),coords(:,2),'k.','markersize',10);

for i=1:length(node)
  nid = i;
  xI  = node(nid,:);
  nI  = normals(nid,:);
  % if norm(nI) > tol, nI = nI / norm(nI); end
  le = 0.01;
  quiver(xI(1), xI(2), gcurvature(nid)*normals(nid,1)/100, ...
                       gcurvature(nid)*normals(nid,2)/100,4,'LineWidth',2,...
         'Color','r','MaxHeadSize',2);  
end

disp([num2str(toc),'   DONE '])
