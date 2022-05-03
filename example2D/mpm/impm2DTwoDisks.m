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
% Two elastic disks come into contact.
% Example taken from Sulsky et al. 1994 paper.
%
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 9 October 2015.
% Status: not working properly for collision (12/10/2015)

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
E      = 1000;              % Young modulus
nu     = 0.3;               % Poisson ratio
rho    = 1000;              % density
v0     = 0.1;
stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);

interval     = 100;% time interval for saving vtp files.
vtkFileName  = 'impm2DTwoDisks';
nodel     = square_node_array([-1 -1],[-1 1],[1 1],[1 -1],3,3);
elemType  = 'Q4';
tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell=0;
lx      = 1;
ly      = 1;
noX0    = 20;      % number of elements along X direction
noY0    = noX0;        % number of elements along Y direction
[mesh]  = buildGrid2D(lx,ly,noX0,noY0, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;
omegaC    = mesh.deltax*mesh.deltay;

%% generate material points
ppc           = [2 2];
cir1.center   = [0.2 0.2];
cir1.radius   = 0.2;
cir2.center   = [0.8 0.8];
cir2.radius   = 0.2;
[res1]        = generateMPForCircle(cir1,ppc,mesh);
[res2]        = generateMPForCircle(cir2,ppc,mesh);

pCount  = 2*size(res1.position,1);
volume  = [res1.volume;res2.volume];
volume0 = volume;
mass    = volume*rho;
coords  = [res1.position;res2.position];
coords0 = coords;
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stress
strain  = zeros(pCount,3);                % strain
velo    = zeros(pCount,2);                % velocity
density = rho*ones(pCount,1);             % density

% initial velocities, initial stress=0
bParticles = [];
for p=1:pCount
  if coords(p,1) < 0.5
    velo(p,:) = [v0 v0];
    di = norm(coords(p,:) - cir1.center);
    if abs(di-0.2)<2e-2
      bParticles = [bParticles;p];
    end
  else
    velo(p,:) = [-v0 -v0];
        di = norm(coords(p,:) - cir2.center);
    if abs(di-0.2)<2e-2
      bParticles = [bParticles;p];
    end
  end
end


%% MLS weight functions
shape = 'circle' ;         % shape of domain of influence
dmax1 = 2.0 ;              % radius = dmax * nodal spacing
dmax2 = 2.0 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;             % using cubic spline weight function

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
plot(coords(:,1),coords(:,2),'k.','markersize',10);
plot(coords(bParticles,1),coords(bParticles,2),'r.','markersize',20);
%plot_mesh(node,element(activeElems,:),'Q4','r-',2.);
plot_mesh(node,element(intersect(find(cells==2),activeElems),:),'Q4','b--',2.3);
axis off

%% nodal quantities
nmass       = zeros(nodeCount,1);  % nodal mass vector
nvelo       = zeros(nodeCount,2);  % nodal velocity vector (final)
nvelo0      = zeros(nodeCount,2);  % nodal velocity vector (begin)
niforce     = zeros(nodeCount,2);  % nodal internal force vector
neforce     = zeros(nodeCount,2);  % nodal external force vector

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.001;
time  = 2;
t     = 0;

istep = 0;
nsteps = floor(time/dtime);
ka     = zeros(nsteps,1);
sa     = zeros(nsteps,1);
ta     = 0:dtime:time;
%%
while ( t < time )
  disp(['time step ',num2str(t)])
  % reset grid data
  nmass(:)     = 0;
  niforce(:)   = 0;
  neforce(:)   = 0;
  nvelo0(:)    = 0;
  % loop over computational cells or elements
  for ie=1:length(activeElems)
    e        = activeElems(ie);
    cellType = cells(e);
    esctr    = element(e,:);
    enode    = node(esctr,:);     % element node coords
    mpts     = mpoints{e};   
    if cellType == 2
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
  nforce    = niforce;
  acce(:,1) = nforce(:,1)./nmass;
  acce(:,2) = nforce(:,2)./nmass;
  if (istep==0), acce = 0.5*acce; end
  nvelo     = nvelo0 + acce*dtime;
  % boundary conditions
  
  
  % update particle velocity, position, stresses
  k = 0; u = 0;
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
      
      dEps           = dtime * 0.5 * (Lp+Lp');
      dsigma         = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
      stress(pid,:)  = stress(pid,:) + dsigma';
      strain(pid,:)  = strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
      
      k = k + 0.5*(velo(pid,1)^2+velo(pid,2)^2)*mass(pid);
      u = u + 0.5*volume(pid)*stress(pid,:)*strain(pid,:)';            
    end
  end

  
  % update the element particle list  
  cells(:) = 0;
  for p=1:pCount
    x = coords(p,1);
    y = coords(p,2);
    e = floor(x/mesh.deltax) + 1 + mesh.numx*floor(y/mesh.deltay);
    pElems(p) = e;
  end
  
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
  
  ka(istep) = k; sa(istep) = u;
end
%% post processing

disp([num2str(toc),'   POST-PROCESSING '])


Ux= zeros(size(node,1),1);
Uy= zeros(size(node,1),1);
sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

VTKPostProcess(node,element,2,'Quad4','../../results/impm/twoDisksGrid',...
    [sigmaXX sigmaYY sigmaXY],[Ux Uy]);
  
mpm=load('mpm-explicit-2disks-dt001.mat');

figure 
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(mpm.ta(1:end),mpm.ka(1:end),'b--','LineWidth',1.6);
plot(ta(1:end),sa(1:end),'r-','LineWidth',2);
plot(mpm.ta(1:end),mpm.sa(1:end),'r--','LineWidth',2);
%plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy (x1E-3)')
legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%axis([0 3 0 3])

%% plot mesh, particles
figure
hold on
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(bodies{1}.elements,:),'Q4','cy-',2.1);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','r-',2.1);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
axis off

disp([num2str(toc),'   DONE '])
