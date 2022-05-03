% This file implements the Material Point Method with CPDI2 interpolation
% described in the article
%
% A. Sadeghirad, R. M. Brannon, and J.E. Guilkey. Second-order convected
% particle domain in- terpolation (CPDI2) with enrichment for weak discontinuities
% at material interfaces. IJNME, 95(11):928-952, 2013.
%
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
%
% Generalized vortex problem.
%
% Vinh Phu Nguyen
% January 2016, Saigon, Vietnam.

%%
addpath ../../fem_util/
addpath ../../fem-functions/
addpath ../../postProcessing/
addpath ../../constitutiveModels/
addpath ../../grid/
addpath ../../util/
addpath ../../basis/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%

E   = 1000;    % Young's modulus
nu  = 0.3;     % Poisson ratio
rho = 1000;    % density
K   = E/3/(1-2*nu);    % bulk modulus
mu    = E/2/(1+nu);% shear modulus
lambda = K - 2/3*mu;
A     = 1/3;

I  = [1 0;0 1];

% be careful with vtkFileName1 and change it according to your computer!!!
interval     = 1;
vtkFileName  = 'vortex';
vtkFileName1 = '../results/cpdi/vortexGrid';
vtkFileName2 = 'vortexExact';


stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])

%%   particle distribution
%

meshFile = 'ring.msh';
pmesh     = load_gmsh (meshFile);

elemType = 'Q4';
numnode  = pmesh.nbNod;
numelem  = pmesh.nbQuads;
node1    = pmesh.POS(:,1:2);
element1 = pmesh.QUADS(1:numelem,1:4);

% store the particle mesh into a structure for convenience
particles.node     = node1;
particles.elem     = element1;
particles.elemType = elemType;

particlesE       = particles;
particlesE.node0 = node1;

pCount  = numelem;                       % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
coords  = zeros(pCount,2);
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stressx
strain  = zeros(pCount,3);                % strain
velo    = zeros(pCount,2);
color   = ones(pCount,1);                 % color, cannot avoid it
u       = zeros(pCount,2);               % displacements of particle corners (visua)

for e = 1:numelem
  coord = node1(element1(e,:),:);
  a     = 0.5*( coord(1,1)*coord(2,2)  - coord(2,1)*coord(1,2) ...
    + coord(2,1)*coord(3,2)  - coord(3,1)*coord(2,2) ...
    + coord(3,1)*coord(4,2)  - coord(4,1)*coord(3,2) ...
    + coord(4,1)*coord(1,2)  - coord(1,1)*coord(4,2) );
  volume(e)  = a;
  mass(e)    = a*rho;
  coords(e,:) = mean(coord); % center of each element=particle
end

volume0 = volume;
coords0 = coords;

%% Computational grid

ghostCell=0;
lx     = 2.7;
ly     = lx;
numx2  = 40;      % number of elements along X direction
numy2  = 40;      % number of elements along Y direction
[mesh]= buildGrid2D(lx,ly,numx2,numy2, ghostCell);

mesh.node(:,1) = mesh.node(:,1) - lx/2;
mesh.node(:,2) = mesh.node(:,2) - lx/2;
element= mesh.element;
node   = mesh.node;

% find boundary nodes

%fixNodes=find(abs(node(:,2)-6*l/2)<1e-10);

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum = zeros(mesh.nodeCount,2);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,2);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,2);  % nodal external force vector
nforce    = zeros(mesh.nodeCount,2);  % nodal force vector

%% plot mesh, particles

hold on
plot_mesh(particles.node,particles.elem,elemType,'black-',1.9);
plot_mesh(node,element,'Q4','r-',0.8);
plot(coords(:,1),coords(:,2),'k.','markersize',20);
%plot(node(fixNodes,1),node(fixNodes,2),'r*','markersize',14);
%axis on
axis off
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'cpdi-ring','-dpdf','-r1000');
%% 

figure(2)
hold on
plot_mesh(node,element,'Q4','r-',0.8);
plot(coords(:,1),coords(:,2),'k.','markersize',14);
lxp = 0.08;
for p=1:numelem
    xc = coords(p,1);
    yc = coords(p,2);
    PG = polyshape([xc-0.5*lxp xc+0.5*lxp xc+0.5*lxp xc-0.5*lxp], [yc-0.5*lxp yc-0.5*lxp yc+0.5*lxp yc+0.5*lxp]);
    plot(PG,'FaceColor','white','EdgeColor', 'black', 'LineWidth', 1.8)
end
axis off
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'gimp-ring','-dpdf','-r1000');


ta = 0;           % time
ka = 0;

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.05*(mesh.deltax/c);
time  = 1/2;
t     = 0;
nsteps = floor(time/dtime);

istep = 1;

while ( t < time )
  disp(['time step ',num2str(t)])
  % reset grid data
  nmass(:)     = 0;
  nmomentum(:) = 0;
  niforce(:)   = 0;
  neforce(:)   = 0;
  % loop over particles
  for p=1:pCount
    sig    = stress(p,:);
    xp     = coords0(p,:);
    % particle mass and momentum to node
    data  = cpdi22D(p,particles,mesh);
    esctr = data.node;
    [bx,by] = vortexBodyForces(xp(1),xp(2),t,mu,rho,A);
    for i=1:length(esctr)
      id              = esctr(i);
      nmass(id)       = nmass(id)       + data.phi(i)*mass(p);
      nmomentum(id,:) = nmomentum(id,:) + data.phi(i)*mass(p)*velo(p,:);
      niforce(id,1)   = niforce(id,1) - volume(p)*(sig(1)*data.dphi(i,1) + sig(3)*data.dphi(i,2));
      niforce(id,2)   = niforce(id,2) - volume(p)*(sig(3)*data.dphi(i,1) + sig(2)*data.dphi(i,2));
      neforce(id,:)   = neforce(id,:) + mass(p)*data.phi(i)*[bx, by];
    end
  end
  
  % update nodal momenta
  nforce    = niforce + neforce;
  %     nforce   (fixNodes,2)  = 0;
  %     nmomentum(fixNodes,2)  = 0;
  
  nmomentum = nmomentum + nforce*dtime;
  
  % update particle velocity and position and stresses
  
  % loop over particles
  for p=1:pCount
    Lp    = zeros(2,2);
    data  = cpdi22D(p,particles,mesh);
    esctr = data.node;
    for i=1:length(esctr)
      id = esctr(i);
      vI = [0 0];
      if nmass(id) > tol
        velo(p,:)  = velo(p,:) + dtime * data.phi(i)*nforce(id,:)/nmass(id);
        vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
      end
      Lp = Lp + vI'*data.dphi(i,:);         % particle gradient velocity
    end
    
    F          = ([1 0;0 1] + Lp*dtime)*reshape(deform(p,:),2,2);
    deform(p,:)= reshape(F,1,4);
    volume(p)  = det(F)*volume0(p);
    J       = det(F);
    b       = F*F';
    sigma   = 1/J*( mu*(b-I) + lambda*log(J)*I );
    stress(p,:)  = [sigma(1,1) sigma(2,2)  sigma(1,2) ];
  end
  
  for c=1:size(particles.node,1)
    xc    = particles.node(c,:);
    ec    = point2ElemIndex(xc,mesh);
    esctr = element(ec,:);
    for i=1:length(esctr)
      id      = esctr(i);
      x       = xc - node(id,:);
      [N,dNdx]= getMPM2D(x,mesh.deltax,mesh.deltay);
      if nmass(id) > tol
        xc = xc + dtime*N*nmomentum(id,:)/nmass(id);
        u(c,:) = dtime*N*nmomentum(id,:)/nmass(id);
      end
    end
    particles.node(c,:) = xc;
  end
  
  for c=1:size(particlesE.node,1)
    xc    = particlesE.node0(c,:);
    R     = sqrt(xc(1)^2+xc(2)^2);
    % (15-32*R+16*R^2)^2 equlas (1-32*(R-1)^2+256*(R-1)^4)!!!
    %alpha = A*0.5*(1-cos(2*pi*t))*(15-32*R+16*R^2)^2;
    alpha = A*0.5*(1-cos(2*pi*t))*(1-32*(R-1)^2+256*(R-1)^4);
    Q = [cos(alpha) -sin(alpha) 0;...
         sin(alpha)  cos(alpha) 0;...
         0            0          1];
    
    x = Q*[xc(1);xc(2);0];
    particlesE.node(c,:) = x(1:2);
  end
  
  % update centroids or particle positions
  for e = 1:numelem
    coord = particles.node(element1(e,:),:);
    a     = 0.5*( coord(1,1)*coord(2,2)  - coord(2,1)*coord(1,2) ...
      + coord(2,1)*coord(3,2)  - coord(3,1)*coord(2,2) ...
      + coord(3,1)*coord(4,2)  - coord(4,1)*coord(3,2) ...
      + coord(4,1)*coord(1,2)  - coord(1,1)*coord(4,2) );
    coords(e,:) = mean(coord); % center of each element=particle
    u(e,:) = coords(e,:) - coords0(e,:);
  end
  
  % VTK output
  
  if (  mod(istep,interval) == 0 )
    vtkFile  = sprintf('../../results/cpdi/%s%d',vtkFileName,istep);
    vtkFile2 = sprintf('../../results/cpdi/%s%d',vtkFileName2,istep);
    data.stress  = [stress zeros(pCount,1)];
    data.pstrain = [];
    data.color   = color;
    data.velo    = velo;
    data.disp    = u;
    VTKParticlesCPDI(particles,vtkFile,data);
    VTKParticlesCPDI(particlesE,vtkFile2,data);
  end
  
  
  % advance to the next time step
  t     = t + dtime;
  istep = istep + 1;
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])


%%
Ux= zeros(size(mesh.node,1),1);
Uy= zeros(size(mesh.node,1),1);
sigmaXX = zeros(size(mesh.node,1),1);
sigmaYY = zeros(size(mesh.node,1),1);
sigmaXY = zeros(size(mesh.node,1),1);

VTKPostProcess(mesh.node,mesh.element,2,'Quad4',vtkFileName1,...
  [sigmaXX sigmaYY sigmaXY],[Ux Uy]);

%%

figure
hold on
plot_mesh(particles.node,particles.elem,elemType,'r-',0.2);
plot_mesh(node,element,'Q4','k-',1.)
plot(coords(:,1),coords(:,2),'k.','markersize',10);
axis off

disp([num2str(toc),'   DONE '])
