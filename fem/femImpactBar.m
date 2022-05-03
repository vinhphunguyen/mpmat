% This file implements the linear FEM.
% The grid consists of four-noded quadrilateral elements.
% Newmark time integration.
%
% Wave propagation in a bar.
%
% Vinh Phu Nguyen
% 1 January 2016.

%%

addpath ../grid/
addpath ../basis/
addpath ../particleGen/
addpath ../constitutiveModels/
addpath ../util/
addpath ../geoMesh/
addpath ../externals/
addpath ../postProcessing/


%%
clc
clear

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

tic;
%%
%
E      = 1;               % Young modulus
nu     = 0.;               % Poisson ratio
rho    = 1;              % density
c      = sqrt(E/rho);
sigma0 = 0.01;

elemType  = 'Q4';
identity  = [1 0;0 1]; 

stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C           = elasticityMatrix(E,nu,stressState);
material.D  = C;

% Newmark time integration constants
beta   = 0.25;
gamma  = 0.5;
CFL    = 0.1;
%% Computational grid
l = 10.;
w = 1.;

noX0      = 200;    % number of elements along X direction
noY0      = 1;      % number of elements along Y direction
ghostCell = 0;
[mesh]    = buildGrid2D(l,w,noX0,noY0, ghostCell);
nodes     = mesh.node;
elements  = mesh.element;
elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
nodes0    = nodes;
mesh.nnode= size(elements,2); % number of node per element
mesh.ndof = size(elements,2)*2; % number of dof per element
mesh.elemType = 'Q4';

%% Boundary nodes (for essential BCs)

fixedNodes  = mesh.rNodes; % nodes
forcedNodes = mesh.lNodes; % nodes
uFixed     = zeros(1,length(fixedNodes));
vFixed     = zeros(1,length(fixedNodes));
udofs      = 2*fixedNodes-1;
vdofs      = 2*fixedNodes;

%%
% Integration points:
noGPs   = 4;
[W,Q]   = quadrature( 2, 'GAUSS', 2 ); % 2x2 Gauss point rule
pCount  = noGPs*elemCount;             % total of Gauss points of the mesh
stress  = zeros(pCount,2,2);           % stress at all GPs, 2x2 matrix for 1 GP
indices = zeros(elemCount,noGPs);
mesh.W  = W;
mesh.Q  = Q;

i = 1;
for e=1:elemCount
  for p=1:length(W)
    indices(e,p) = i; i = i + 1;
  end
end

figure('Color','white')
hold on
plot_mesh(nodes,elements,elemType,'k-',1.);
plot(nodes(fixedNodes,1),nodes(fixedNodes,2),'bs','MarkerSize',9,'MarkerFaceColor','blue')
axis off

%% nodal quantities, matrices...
ndof      = 2;  
dofCount  = nodeCount*ndof;
mesh.dofCount = dofCount;

% nodal accelerations, velocities, displacements at time 't'
nacce0    = zeros(dofCount,1);  % nodal  acceleration
nvelo0    = zeros(dofCount,1);  % nodal  velocity vector 
ndisp0    = zeros(dofCount,1);  % nodal displacement vector
% nodal accelerations, velocities, displacements at time 't+dtime'
nacce     = zeros(dofCount,1);  % nodal  acceleration
nvelo     = zeros(dofCount,1);  % nodal  velocity vector 
ndisp     = zeros(dofCount,1);  % nodal displacement vector

massMat   = zeros(dofCount,dofCount); % consistent mass matrix
fext      = zeros(dofCount,1);        % external force vector

sigma     = zeros(nodeCount,3);  % nodal external force vector

%Compute consistent mass matrix once
for e=1:elemCount
  esctr = elements(e,:);
  enode = nodes(esctr,:);
  nn    = length(esctr);
  sctr(1,1:2:2*nn) = 2*esctr-1;    
  sctr(1,2:2:2*nn) = 2*esctr  ;
  % loop over Gauss point
  for p=1:length(W)
    pt       = Q(p,:);
    [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
    J0       = enode'*dNdxi;                  % element Jacobian matrix
    detJ     = det(J0);
    Nmat(1,1:2:2*nn)   = N;
    Nmat(2,2:2:2*nn)   = N;    
    massMat(sctr,sctr) = massMat(sctr,sctr) + Nmat' * Nmat * rho * detJ * W(p);
  end
end

%% Time loop
dtime  = 0.01;%CFL*mesh.deltax/c;
dtime2 = dtime*dtime;
time   = 500*dtime;0.5*l/c;
t      = 0.;
istep  = 0;
interval = 10;

% initialise nodal velocities and displacements
% for i=1:nodeCount
%   nvelo(i,1) = mmsV1(nodes(i,1),0);
%   nvelo(i,2) = mmsV2(nodes(i,2),0);
%   
%   ndisp(i,1) = mmsU1(nodes(i,1),0);
%   ndisp(i,2) = mmsU2(nodes(i,2),0);
% end

nsteps = floor(time/dtime);
pDisp  = zeros(nsteps,1);
ta     = 0:dtime:time;
tol    = 1e-5;

markedElem = noX0;

while ( t < time )
  dtilde = ndisp + dtime*nvelo + 0.5*dtime2*(1-2*beta)*nacce;   
  % compute the stiffness matrix
  [stiffMat] = computeStiffnessMatrix(mesh,material);
  % compute the external force
  fext(2*forcedNodes-1) = 0.5*sigma0;
  % modified stiffness matrix
  K      = 1/(beta*dtime2)*massMat + stiffMat;
  % Residual vector
  res    = 1/(beta*dtime2)*massMat*dtilde + fext;
  % Boundary conditions
  [K,res]  = applyDirichletBCs(K,res,udofs,vdofs,uFixed,vFixed);
  % Solving for displacement increment
  ndisp    = K\res;
  % update accelerations/velocities 
  nacce    = 1/(beta*dtime2)*(ndisp-dtilde);
  nvelo    = nvelo0 + (1-gamma)*dtime*nacce0 + gamma*dtime*nacce;
  
  % advance to the next time step
  t        = t + dtime;
  istep    = istep + 1;
  
  nacce0   = nacce;
  nvelo0   = nvelo;
  ndisp0   = ndisp;
  
  
%   if (  mod(istep-1,interval) == 0 )   
%     vtuFile = sprintf('%s%d','femImpVibratingCantilever',istep-1);
%     VTKPostProcess(nodes,elements,2,'Quad4',vtuFile,sigma,...
%       [ndisp(1:2:dofCount) ndisp(2:2:dofCount)]);
%   end
end
%%
disp([num2str(toc),'   DONE ']);

numSig = zeros(1*mesh.elemCount,1);
numLoc = zeros(1*mesh.elemCount,1);
ii     = 1;
for e=1:mesh.elemCount
  esctr = mesh.element(e,:);
  enode = mesh.node(esctr,:);
  nn    = length(esctr);
  sctr(1,1:2:2*nn) = 2*esctr-1;    
  sctr(1,2:2:2*nn) = 2*esctr  ;  
  eDisp            = ndisp(sctr);
  % loop over Gauss point
  for p=1:1
    pt       = mesh.Q(p,:);
    [N,dNdxi]= lagrange_basis(mesh.elemType,pt);   % element shape functions
    J0       = enode'*dNdxi;                       % element Jacobian matrix
    detJ     = det(J0);
    dNdx     = dNdxi/J0;  
    B(1,1:2:2*nn)  = dNdx(:,1)';
    B(2,2:2:2*nn)  = dNdx(:,2)';
    B(3,1:2:2*nn)  = dNdx(:,2)';
    B(3,2:2:2*nn)  = dNdx(:,1)';    
    sig            = C*B*eDisp;    xx = N'*enode;
    numSig(ii)     = sig(1);
    numLoc(ii)     = xx(1);
    ii             = ii + 1;
  end
end

% exact solution
x        = linspace(0,l,100);
exactSig = zeros(length(x),1);
for i=1:length(x)
  exactSig(i) = -sigma0*heaviside(c*time-x(i));
end

%%
figure
hold on
plot(numLoc,numSig,'red-','LineWidth',1.8)
plot(x,exactSig,'blue-','LineWidth',1.8)
xlabel('x/L')
ylabel('stress')
set(gca,'FontSize',16)
legend('FEM','Exact')
grid on

%%
