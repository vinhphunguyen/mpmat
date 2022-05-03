% This file implements the Improved Material Point Method.
% The grid: four-noded bilinear elements.
% Leapfrog time integration.
%
% Moving Least Square approximation used to construct particle data on the
% grid (nodes/centers).
% MLS introduces two parameters into the problem: smoothing length (domain of influence).
% Generally, vvelocity MLS has one smoothing length and density/stress MLS have
% another smoothing length. Parameter study is needed.
%
% Large deformation of a vibrating cantilever beam.
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 9 October 2015.

%%
%parpool('open',8);

%%

addpath ../../grid/
addpath ../../basis/
addpath ../../particleGen/
addpath ../../constitutiveModels/
addpath ../../util/
addpath ../../geoMesh/
addpath ../../externals/
addpath ../../postProcessing/
addpath ../../mls/

%%
clc
clear
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
E      = 1e6;               % Young modulus
nu     = 0.3;               % Poisson ratio
rho    = 1050;              % density
lambda = E*nu/(1+nu)/(1-2*nu);
mu     = E/2/(1+nu);
c      = sqrt(E/rho);
g      = -10;

identity  = [1 0;0 1]; 

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell=0;
lx      = 6;
ly      = 8;
noX0    = 24;      % number of elements along X direction
noY0    = 32;        % number of elements along Y direction
[mesh]  = buildGrid2D(lx,ly,noX0,noY0, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;
omegaC    = mesh.deltax*mesh.deltay;

%% generate material points
ppc           = [3 3];
square.x      = [0 4];
square.y      = [0 1];
[pmesh]       = buildGrid2D(4,1,12,3, ghostCell);
[res]         = generateMPForRectangle(square,ppc,pmesh);
res.position(:,2) = res.position(:,2) + 4;

pCount  = size(res.position,1);
volume  = res.volume;
volume0 = res.volume;
mass    = res.volume*rho;
coords  = res.position;
coords0 = coords;
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stress
velo    = zeros(pCount,2);                % velocity
density = rho*ones(pCount,1);             % density

% initial velocities, initial stress=0
% for p=1:pCount
%   velo(p,1) = ...;
%   velo(p,2) = ...;
% end

idx = intersect(find(abs(coords(:,1)-4)<mesh.deltax),...
                find(abs(coords(:,2)-4)<mesh.deltay));
              
marked = 134;              

%% MLS weight functions
shape = 'circle' ;         % shape of domain of influence
dmax1 = 3.5 ;              % radius = dmax * nodal spacing
dmax2 = 3.1 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;             % using cubic spline weight function

% Domain of influence for every nodes(sampling points)
% Uniformly distributed nodes
% Definition : rad = dmax*deltaX
di1 = ones(length(coords),1)*dmax1*mesh.deltax;
di2 = ones(length(coords),1)*dmax2*mesh.deltax;

%% plot mesh, particles
figure(1)
hold on
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(bodies{1}.elements,:),'Q4','cy-',2.1);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','r-',2.1);
plot(coords(:,1),coords(:,2),'r.','markersize',40);
plot(coords(46,1),coords(46,2),'rs','markersize',40);
axis off

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

pElems  = ones(pCount,1);
mpoints = cell(mesh.elemCount,1);

for p=1:pCount
  x = coords(p,1);
  y = coords(p,2);
  e = floor(x/mesh.deltax) + 1 + mesh.numx*floor(y/mesh.deltay);
  pElems(p) = e;
end

for e=1:mesh.elemCount
  id  = find(pElems==e);
  mpoints{e}=id;
end


activeElems = unique(pElems);
activeNodes = unique(element(activeElems,:));

%% nodal quantities
nmass       = zeros(nodeCount,1);  % nodal mass vector
nvelo       = zeros(nodeCount,2);  % nodal velocity vector (final)
nvelo0      = zeros(nodeCount,2);  % nodal velocity vector (begin)
niforce     = zeros(nodeCount,2);  % nodal internal force vector
neforce     = zeros(nodeCount,2);  % nodal external force vector

cellDensity = zeros(elemCount,1);  % cell-centerd density
cellStress  = zeros(elemCount,3);  % cell-centerd stress

for ie=1:length(activeElems)
  e     = activeElems(ie);
  esctr = element(e,:);
  enode = node(esctr,:);
  xc    = mean(enode);
  index = defineSupport2D(coords,xc,di2);
  %if length(index) <= 2, disp('A singular'); end
  phi   = mlsLinearBasis2D(xc,index,coords,di2,form);
  cellDensity(e) = cellDensity(e) + dot(phi,density(index));
  %cellStress(e)  = cellStress(e)  + dot(phi,s(index));
end

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.1*mesh.deltax/c;%0.002;
time  = 3;
t     = 0;

istep = 1;
nsteps = floor(time/dtime);
pDisp  = zeros(nsteps,1);
ta     = 0:dtime:time;

while ( t < time )
  disp(['time step ',num2str(t)])
  % reset grid data
  nmass(:)     = 0;
  niforce(:)   = 0;
  neforce(:)   = 0;
  nvelo0(:)    = 0;
  % loop over computational cells or elements
  for ie=1:length(activeElems)
    e     = activeElems(ie);    
    esctr = element(e,:);
    enode = node(esctr,:);     % element node coords
    mpts  = mpoints{e};
    % nodal mass and internal force with 1 point quadrature
    xc    = mean(enode);
    sigc  = cellStress(e,:);
    rhoc  = cellDensity(e);
    for i=1:length(esctr)
      id      = esctr(i);
      x       = xc - node(id,:);
      [N,dNdx]= getMPM2D(x,mesh.deltax,mesh.deltay);
      %nmass(id)       = nmass(id)       + N*rhoc*omegaC;
      niforce(id,1)   = niforce(id,1) - omegaC*(sigc(1)*dNdx(1) + sigc(3)*dNdx(2));
      niforce(id,2)   = niforce(id,2) - omegaC*(sigc(3)*dNdx(1) + sigc(2)*dNdx(2));
      %neforce(id,2)   = neforce(id,2) + omegaC*rhoc*N*g;
    end
    % external force and mass
    for p=1:length(mpts)
        pid  = mpts(p);
        xx   = coords(pid,:);
        XX   = coords0(pid,:);
        mm   = mass(pid);
        for i=1:length(esctr)
            id              = esctr(i);
            x               = xx - node(id,:);
            [N,dNdx]        = getMPM2D(x,mesh.deltax,mesh.deltay);  
            nmass(id)       = nmass(id)     + N*mm;
            neforce(id,2)   = neforce(id,2) + mm*N*g;
        end
    end
    
  end
  
  % project velocity to grid nodes (MLS)
  for ii=1:length(activeNodes)
    i     = activeNodes(ii);
    pt    = node(i,:);
    index = defineSupport2D(coords,pt,di1);
    %if length(index) <= 2, disp('A singular'); end
    phi   = mlsLinearBasis2D(pt,index,coords,di1,form);
    nvelo0(i,1) = nvelo0(i,1) + dot(phi,velo(index,1));
    nvelo0(i,2) = nvelo0(i,2) + dot(phi,velo(index,2));
  end
  
  % update nodal velocity
  nforce    = niforce + neforce;
  acce(:,1) = nforce(:,1)./nmass;
  acce(:,2) = nforce(:,2)./nmass;
  if (istep==0), acce = 0.5*acce; end
  nvelo     = nvelo0 + acce*dtime;
  % boundary conditions
  nvelo(mesh.lNodes,:)  = 0.;
  
  % update particle velocity, position, stresses
  for ie=1:length(activeElems)
    e     = activeElems(ie);   
    esctr = element(e,:);
    enode = node(esctr,:);
    mpts  = mpoints{e};   
    % loop over particles
    for p=1:length(mpts)
      pid  = mpts(p);
      xp   = coords(pid,:);
      Lp   = zeros(2,2);
      for i=1:length(esctr)
        id = esctr(i);
        vI = [0 0];
        x       = xp - node(id,:);
        [N,dNdx]= getMPM2D(x,mesh.deltax,mesh.deltay);
        if nmass(id) > 0
          velo(pid,:)   = velo(pid,:)   + N*(nvelo(id,:)-nvelo0(id,:));
          coords(pid,:) = coords(pid,:) + dtime*N*nvelo(id,:);
          vI            = nvelo(id,:);  % nodal velocity  
        end
        Lp = Lp + vI'*dNdx;         % particle gradient velocity
      end
      
      F             = ([1 0;0 1] + Lp*dtime)*reshape(deform(pid,:),2,2);
      deform(pid,:) = reshape(F,1,4);
      volume(pid)   = det(F)*volume0(pid);
      J             = det(F);
      if (J<0), disp('error');end;
      density(pid)  = rho/J;      
      b             = F*F';
      sigma         = 1/J*( mu*(b-identity) + lambda*log(J)*identity );
      stress(pid,:)  = [sigma(1,1) sigma(2,2)  sigma(1,2)];
    end
  end

  
  % update the element particle list  
  for p=1:pCount
    x = coords(p,1);
    y = coords(p,2);
    e = floor(x/mesh.deltax) + 1 + mesh.numx*floor(y/mesh.deltay);
    pElems(p) = e;
  end
  
  for e=1:mesh.elemCount
    id  = find(pElems==e);
    mpoints{e}=id;
  end
  
  activeElems = unique(pElems);
  activeNodes = unique(element(activeElems,:));
      
  % project particle density/stress to grid centers (MLS)
  cellDensity(:) = 0;
  cellStress(:)  = 0;
  for ie=1:length(activeElems)
    e     = activeElems(ie);   
    esctr = element(e,:);
    enode = node(esctr,:);        
    xc    = mean(enode);    
    index = defineSupport2D(coords,xc,di2);
    %if length(index) <= 2, disp('A singular'); end
    phi   = mlsLinearBasis2D(xc,index,coords,di2,form);
    cellDensity(e)   = cellDensity(e)   + dot(phi,density(index));
    cellStress(e,1)  = cellStress(e,1)  + dot(phi,stress(index,1));
    cellStress(e,2)  = cellStress(e,2)  + dot(phi,stress(index,2));
    cellStress(e,3)  = cellStress(e,3)  + dot(phi,stress(index,3));
  end
  
  t     = t + dtime;
  istep = istep + 1;  
  
  pDisp(istep) = coords(marked,2)-coords0(marked,2);
end
%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

figure
hold on
plot(ta(2:end),pDisp(1:end),'blue-','LineWidth',1.8)
xlabel('time')
ylabel('displacement')
set(gca,'FontSize',16)
%legend('G=0.0001','G=0.01', 'G=0.05')
grid on

disp([num2str(toc),'   DONE '])
