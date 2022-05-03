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
% Not working: numerical fracture near support!!!
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

elemType  = 'Q4';
identity  = [1 0;0 1]; 
nodel     = square_node_array([-1 -1],[-1 1],[1 1],[1 -1],3,3);

interval     = 10;% time interval for saving vtp files.
vtkFileName  = 'impm2DVibratingCantilever';

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell=0;
lx      = 6;
ly      = 8;
noX0    = 36;      % number of elements along X direction
noY0    = 48;        % number of elements along Y direction
[mesh]  = buildGrid2D(lx,ly,noX0,noY0, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;
omegaC    = mesh.deltax*mesh.deltay;
nodalSup  = mesh.nodalSup;

%% generate material points
ppc           = 3;
square.x      = [0 4];
square.y      = [0 1];
%[pmesh]       = buildGrid2D(4,1,12,3, ghostCell);
%[res]         = generateMPForRectangle(square,ppc,pmesh);
%res.position(:,2) = res.position(:,2) + 4;
% particles with ones on the boundaries (different from most MPM)
[pmesh]         = buildGrid2D(4,1,32,8, ghostCell);
pmesh.node(:,2) = pmesh.node(:,2) + 4;

% no longer needs particle volumes/mass
pCount  = size(pmesh.node,1);
%volume  = res.volume;
%volume0 = res.volume;
%mass    = res.volume*rho;
coords  = pmesh.node;
coords0 = coords;
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stress
velo    = zeros(pCount,2);                % velocity
density = rho*ones(pCount,1);             % density
bParticles = unique([pmesh.bNodes;pmesh.rNodes;pmesh.tNodes]);

% initial velocities, initial stress=0
% for p=1:pCount
%   velo(p,1) = ...;
%   velo(p,2) = ...;
% end

idx = intersect(find(abs(coords(:,1)-4)<mesh.deltax),...
                find(abs(coords(:,2)-4)<mesh.deltay));
              
marked = 33;        % number of elements along X (pmesh) + 1

%% MLS weight functions
shape = 'circle' ;         % shape of domain of influence
dmax1 = 3.0 ;              % radius = dmax * nodal spacing
dmax2 = 3.0 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;   % using cubic spline weight function

% Domain of influence for every nodes(sampling points)
% Uniformly distributed nodes
% Definition : rad = dmax*deltaX
di1 = ones(length(coords),1)*dmax1*mesh.deltax;
di2 = ones(length(coords),1)*dmax2*mesh.deltax;

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

cells       = zeros(elemCount,1);  % classify cells
%cells(c)   = 0: void cell
%cells(c)   = 1: full cell
%cells(c)   = 2: cut cell

pElems  = ones(pCount,1);
mpoints = cell(elemCount,1);

for p=1:pCount
  x = coords(p,1);
  y = coords(p,2);
  e = floor(x/mesh.deltax) + 1 + mesh.numx*floor(y/mesh.deltay);
  pElems(p) = e;
end

for e=1:elemCount
  id  = find(pElems==e);
  mpoints{e}=id;
  % detect cut cells
  if ~isempty(intersect(bParticles,id)) 
    cells(e) = 2;
  end
end


activeElems = unique(pElems);
activeNodes = unique(element(activeElems,:));


%% plot mesh, particles
figure(1)
hold on
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(bodies{1}.elements,:),'Q4','cy-',2.1);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','r-',2.1);
plot(coords(:,1),coords(:,2),'b.','markersize',30);
plot(coords(marked,1),coords(marked,2),'cy.','markersize',20);
%plot(coords(46,1),coords(46,2),'bs','markersize',10);
%plot_mesh(node,element(activeElems,:),'Q4','r-',2.);
plot_mesh(node,element(nodalSup{241},:),'Q4','r-',2.);
plot_mesh(node,element(cells==2,:),'Q4','b--',1.3);
%plot(node(activeNodes,1),node(activeNodes,2),'c.','markersize',18);
axis off


%% nodal quantities
nmass       = zeros(nodeCount,1);  % nodal mass vector
nvelo       = zeros(nodeCount,2);  % nodal velocity vector (final)
nvelo0      = zeros(nodeCount,2);  % nodal velocity vector (begin)
niforce     = zeros(nodeCount,2);  % nodal internal force vector
neforce     = zeros(nodeCount,2);  % nodal external force vector

[W,Q]       = quadrature(2,'GAUSS',2);      

%% Solver
disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.002;
time  = 1.6;
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
      
  % project velocity to grid nodes (MLS)
  for ii=1:length(activeNodes)
    i     = activeNodes(ii);
    % if boundary node, then it might have a zero mass (in the sense of standard MPM)
    sup   = nodalSup{i};
    if any(cells(sup)==2)
      mm = 0;
      for ic=1:length(sup)
        cell  = sup(ic);
        mpts  = mpoints{cell};
        % loop over particles
        for p=1:length(mpts)
          pid  = mpts(p);
          x       = coords(pid,:) - node(i,:);
          [N,dNdx]= getMPM2D(x,mesh.deltax,mesh.deltay);
          mm      = mm + N;
        end
      end
      if mm < 1e-6 
        cells(sup) = 1;
        continue;
      end
    end
    pt    = node(i,:);
    index = defineSupport2D(coords,pt,di1);
    %if length(index) <= 2, disp('A singular'); end
    phi   = mlsLinearBasis2D(pt,index,coords,di1,form);
    nvelo0(i,1) = nvelo0(i,1) + dot(phi,velo(index,1));
    nvelo0(i,2) = nvelo0(i,2) + dot(phi,velo(index,2));
  end
  
  % loop over computational cells or elements
  for ie=1:length(activeElems)
    e        = activeElems(ie);
    cellType = cells(e);
    esctr    = element(e,:);
    enode    = node(esctr,:);     % element node coords
    mpts     = mpoints{e};   
    if     cellType == 1
      continue;
    elseif cellType == 2
      [aa] = hierarchicalGaussQuad(2,enode,nodel,mpts,coords,1,1);
      W    = aa.W; Q = aa.Q;
    else
      [W,Q] = quadrature(2,'GAUSS',2);      
    end
    % nodal mass and internal force with full quadrature
    for p=1:length(W)
      pt       = Q(p,:);
      [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
      J0       = enode'*dNdxi;                  % element Jacobian matrix
      dNdx     = dNdxi/J0;                      
      wt       = W(p)*det(J0); 
      xp       = N'*enode;            
      index    = defineSupport2D(coords,xp,di2);
      %if length(index) <= 2, disp('A singular'); end
      phi      = mlsLinearBasis2D(xp,index,coords,di2,form);
      rhoc     = dot(phi,density(index));
      sigc(1)  = dot(phi,stress(index,1));
      sigc(2)  = dot(phi,stress(index,2));
      sigc(3)  = dot(phi,stress(index,3));
      
      nmass(esctr)       = nmass(esctr)     + N*rhoc*wt;
      niforce(esctr,:)   = niforce(esctr,:) - wt*dNdx*[sigc(1) sigc(3);sigc(3) sigc(2)];      
      neforce(esctr,2)   = neforce(esctr,2) + wt*rhoc*N*g;      
    end
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
%      volume(pid)   = det(F)*volume0(pid);
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
  cells(:) = 0;
  for e=1:mesh.elemCount
    id  = find(pElems==e);
    mpoints{e}=id;
    % detect cut cells
    if ~isempty(intersect(bParticles,id))
      cells(e) = 2;
    end
  end
  
  activeElems = unique(pElems);
  activeNodes = unique(element(activeElems,:));
       
  t     = t + dtime;
  istep = istep + 1;  
  
  if (  mod(istep,interval) == 0 )   
    data.stress  = [stress zeros(pCount,1)];
    vtkFile = sprintf('../../results/impm/%s%d',vtkFileName,istep);
    VTKParticles(coords,vtkFile,data);
    %vtkFileName1 = sprintf('../results/%s%d','mpmTwoDisksGrid',istep);
    %VTKPostProcess(node,element,2,'Quad4',vtkFileName1,gStress', gDisp');
  end
  
  pDisp(istep) = coords(marked,2)-coords0(marked,2);
end
%% post processing


Ux= zeros(size(node,1),1);
Uy= zeros(size(node,1),1);
sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

VTKPostProcess(node,element,2,'Quad4','../../results/impm/vibratingCantileverGrid',...
    [sigmaXX sigmaYY sigmaXY],[Ux Uy]);

%%
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
