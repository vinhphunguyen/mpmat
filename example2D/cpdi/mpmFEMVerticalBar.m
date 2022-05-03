% This file implements the  coupled FEM-Material Point Method
%
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
%
% Vertical bar with extreme deformation.
%
% Vinh Phu Nguyen
% Monash University, Australia
% January 2016

%%

addpath ../../fem_util/
addpath ../../fem-functions/
addpath ../../post-processing/
addpath ../../constitutiveModels/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%

E      = 1;           % Young's modulus
nu     = 0.3;         % Poisson ratio
rho    = 1050e-12;    % density
K      = E/3/(1-2*nu);  % bulk modulus
mu     = E/2/(1+nu);    % shear modulus
lambda = K - 2/3*mu;

g     = 700e3; % gravity
bodyf = [0 -g];

I  = [1 0;0 1];

% be careful with vtkFileName1 and change it according to your computer!!!
interval     = 1;
vtkFileName  = 'mpmfem-bar';
vtkFileName1 = '../results/cpdi2/verticalBar/barGrid';

tic;

disp([num2str(toc),'   INITIALISATION '])

%%   particle distribution
%

l     = 1000;
numx2 = 6;      % number of elements along X direction
numy2 = 6;      % number of elements along Y direction
[pmesh]= buildGrid2D(l,l,numx2,numy2, 0);

pmesh.node(:,1) = pmesh.node(:,1) +   l/2;
pmesh.node(:,2) = pmesh.node(:,2) + 4*l/2;
pmesh.elem = pmesh.element;

% store the particle mesh into a structure for convenience
elemType           = 'Q4';
pmesh.elemType     = 'Q4';

pCount  = pmesh.nodeCount;                % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
coords0 = pmesh.node;
coords  = coords0;
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stressx
strain  = zeros(pCount,3);                % strain
velo    = zeros(pCount,2);                % velocities
displ   = zeros(pCount,2);                % displacements
fint    = zeros(pCount,2);                % internal forces (FEM)
centroid= zeros(pmesh.elemCount,2);   % centroids of the particle elements
deformC = repmat([1 0 0 1],pmesh.elemCount,1);     % gradient deformation

[W,Q] = quadrature(2, 'GAUSS', 2); 
%W=4;Q=[0 0];

% determine particle mass
M = zeros(pCount,pCount);
for e=1:pmesh.elemCount
  sctr   = pmesh.element(e,:);         % element scatter vector
  nn     = length(sctr);
  pts    = pmesh.node(sctr,:);          % element nodes
  for gp=1:size(W,1)                    % loop over Gauss points
    pt      = Q(gp,:);
    wt      = W(gp);
    [N,dNdxi]=lagrange_basis(pmesh.elemType,pt);   % element shape functions
    dxdxi = pts'*dNdxi;                  % Jacobian matrix
    detJ  = det(dxdxi);
    mm    = N * N' * rho * detJ * wt;
    M(sctr,sctr) = M(sctr,sctr) + mm;
  end
end
% lumped mass matrix
for p=1:pCount
  mass(p) = sum(M(p,:));
end


%% Computational grid

ghostCell=0;
lx     = (l/2)*4;
ly     = (l/2)*7;
numx2  = 4;      % number of elements along X direction
numy2  = 7;      % number of elements along Y direction
[mesh]= buildGrid2D(lx,ly,numx2,numy2, ghostCell);
element= mesh.element;
node   = mesh.node;

% find boundary nodes

fixNodes=find(abs(node(:,2)-6*l/2)<1e-10);

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

pElems  = ones(pCount,1);
mpoints = cell(mesh.elemCount,1);
xmin = 0; ymin = 0;
for p=1:pCount
  x = coords(p,1);
  y = coords(p,2);
  e = floor((x-xmin)/mesh.deltax) + 1 + numx2*floor((y-ymin)/mesh.deltay);
  pElems(p) = e;
end

for e=1:mesh.elemCount
  id  = find(pElems==e);
  mpoints{e}=id;
end

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum0= zeros(mesh.nodeCount,2);  % nodal momentum vector
nmomentum = zeros(mesh.nodeCount,2);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,2);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,2);  % nodal external force vector
nforce    = zeros(mesh.nodeCount,2);  % nodal force vector

%% plot mesh, particles
figure(1)
hold on
plot_mesh(pmesh.node,pmesh.element,elemType,'black-',1.2);
plot_mesh(mesh.node,mesh.element,'Q4','r-',1.2);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
plot(node(fixNodes,1),node(fixNodes,2),'r*','markersize',14);
%plot_mesh(pmesh.node,pmesh.element(3,:),elemType,'red-',2);
axis off

ta = 0;           % time
ka = 0;

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.1*(mesh.deltax/c);
time  = 0.11;
t     = 0;
nsteps = floor(time/dtime);

istep = 1; 

% one-point quadrature
%W=4;Q=[0 0];

while ( t < time )
  disp(['time step ',num2str(t)])
  % reset grid data
  nmass(:)     = 0;
  nmomentum0(:)= 0;
  niforce(:)   = 0;
  neforce(:)   = 0;
  
  for e=1:mesh.elemCount      % loop over computational cells or elementss
    esctr = mesh.element(e,:);      % element connectivity
    enode = mesh.node(esctr,:);     % element node coords
    mpts  = mpoints{e};        % particles inside element e
    for p=1:length(mpts)       % loop over particles
      pid    = mpts(p);
      xxp    = coords(pid,:);
      Mp     = mass(pid);
      vp     = velo(pid,:);
      forces = fint(pid,:);
      
      pt(1)= (2*xxp(1)-(enode(1,1)+enode(2,1)))/mesh.deltax;
      pt(2)= (2*xxp(2)-(enode(2,2)+enode(3,2)))/mesh.deltay;
      [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions            
      % loop over nodes of current element "ie"
      for i=1:length(esctr)
        id               = esctr(i);
        Ni               = N(i);
        nmass(id)        = nmass(id)         + Ni*Mp;
        nmomentum0(id,:) = nmomentum0(id,:)  + Ni*Mp*vp;
        niforce(id,:)    = niforce(id,:)     - Ni*forces;
        neforce(id,:)    = neforce(id,:)     + Ni*Mp*bodyf;
      end
    end
  end
  
  % update nodal momenta
  nforce    = niforce + neforce;
  nmomentum = nmomentum0 + nforce*dtime;
  nmomentum(fixNodes,2)  = 0;
  
  %% update particle velocity and positions for all particles
  for e=1:mesh.elemCount
    esctr = mesh.element(e,:);      % element connectivity
    enode = mesh.node(esctr,:);     % element node coords
    mpts  = mpoints{e};             % particles inside element e
    % loop over particles
    for p=1:length(mpts)
      pid  = mpts(p);
      xp   = coords(pid,:);
      xp0  = coords0(pid,:);
      Mp   = mass(pid);
      vp   = velo(pid,:);
      up   = displ(pid,:);
      % retrieve the grid functions/grads
      %N    = gBasis(ib,pid,:);
      for i=1:length(esctr)
        id      = esctr(i);
        x       = xp0 - node(id,:);
        [N,~]   = getMPM2D(x,mesh.deltax,mesh.deltay);
        if nmass(id) > 0
          massInv  = (1/nmass(id))*N;
          vI       = nmomentum(id,:)/nmass(id);   % nodal velocity
          vp       = vp +  (nmomentum(id,:)-nmomentum0(id,:))*massInv;
          xp       = xp + dtime * nmomentum(id,:)*massInv;
          up       = up + dtime * vI*N;
        end
      end
      velo(pid,:)   = vp;
      coords(pid,:) = xp;
      displ(pid,:)  = up;
    end
  end
  
  pmesh.node = coords;
  % FEM: update GP stress and then particle internal forces
  f    = zeros(pCount*2,1);
  for e=1:pmesh.elemCount
    sctr   = pmesh.element(e,:);         %  element scatter vector
    nn     = length(sctr);
    sctrB(1,1:nn)      = sctr;
    sctrB(1,nn+1:2*nn) = sctr+pCount;
    B      = zeros(3,2*nn);
    pts    = pmesh.node(sctr,:);          % element nodes' coords
    ue     = displ(sctr,:);   
    ue     = ue';    
    
    for gp=1:size(W,1)                    % loop over Gauss points
      pt      = Q(gp,:);
      wt      = W(gp);
      [~,dNdxi]=lagrange_basis(pmesh.elemType,pt);   % element shape functions
      % Jacobian matrix
      dxdxi = pts'*dNdxi;
      dNdx  = dNdxi/dxdxi;
      detJ  = det(dxdxi);
      % B matrix
      B(1,1:nn)      = dNdx(:,1)';
      B(2,nn+1:2*nn) = dNdx(:,2)';
      B(3,1:nn)      = dNdx(:,2)';
      B(3,nn+1:2*nn) = dNdx(:,1)';
      
      
      F        = inv(I - ue*dNdx);       
      invF     = inv(F);
      detF     = det(F);      
      if detF < 0, error('negative J'); end
      sig      = (mu*(F*F'-I) + lambda*log(detF)*I)/detF;      
      % internal force
      f(sctrB) = f(sctrB) + B'*[sig(1,1);sig(2,2);sig(1,2)]*wt*detJ;           
    end
  end
  fint = reshape(f,pCount,2);
  
  
  % VTK output
  
      if (  mod(istep,interval) == 0 )
          vtkFile = sprintf('../../results/cpdi/%s%d',vtkFileName,istep);
          data.stress  = [stress zeros(pCount,1)];
          data.pstrain = [];
          data.velo    = velo;
          VTKParticlesCPDI(pmesh,vtkFile,data);
      end
  
  % update particle element 
 
  for p=1:pCount
    x = coords(p,1);
    y = coords(p,2);
    e = floor((x-xmin)/mesh.deltax) + 1 + numx2*floor((y-ymin)/mesh.deltay);
    pElems(p) = e;
  end
  
  for e=1:mesh.elemCount
    id  = find(pElems==e);
    mpoints{e}=id;
  end

  % advance to the next time step
  
  t     = t + dtime;
  istep = istep + 1;
  
  % store time,velocty for plotting
  elem = 3;
  sctr = mesh.element(elem,:);
  ue     = pmesh.node(sctr,:)-coords0(sctr,:); 
  ueY    = ue(:,2);
  [N,~]=lagrange_basis(pmesh.elemType,[0 0]);   % element shape functions
  u  = dot(N,ueY);
  ta = [ta;t];
  ka = [ka;u];
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])


ss=load('cpdi2VerticalBar.mat');

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b*','LineWidth',1.6);
plot(ss.ta(1:end),ss.ka(1:end),'r-','LineWidth',1.6);
xlabel('Time')
ylabel('Displacement')
legend('CPDI2','CPDI')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%axis([0 0.045 0 -1800])

%%

figure
hold on
plot_mesh(pmesh.node,pmesh.element,elemType,'black-',1.2);
plot_mesh(mesh.node,mesh.element,'Q4','r-',1.2);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
%plot(node(fixNodes,1),node(fixNodes,2),'r*','markersize',14);
%plot_mesh(pmesh.node,pmesh.element(3,:),elemType,'red-',2);
axis on

disp([num2str(toc),'   DONE '])
