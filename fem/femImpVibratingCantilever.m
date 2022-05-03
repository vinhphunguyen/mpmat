% This file implements the implicit Updated Lagrangian FEM.
% The grid consists of four-noded quadrilateral elements.
% Newmark time integration.
%
% Large deformation of a vibrating cantilever beam.
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
E      = 1e6;               % Young modulus
nu     = 0.3;               % Poisson ratio
rho    = 1050;              % density
lambda = E*nu/(1+nu)/(1-2*nu);
mu     = E/2/(1+nu);
c      = sqrt(E/rho);
g      = -10;

elemType  = 'Q4';
identity  = [1 0;0 1]; 
% Newmark time integration constants
beta   = 0.25;
gamma  = 0.5;

material.mu     = mu;
material.lambda = lambda;
%% Computational grid
l = 4.;
w = 1.;

noX0      = 12;     % number of elements along X direction
noY0      = 3;      % number of elements along Y direction
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

fixedNodes = find(nodes(:,1)<1e-12); % nodes
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
for e=1:elemCount
  esctr  = elements(e,:);
  enode  = nodes(esctr,:);    
  % loop over integration points
  for p=1:length(W)
    pt       = Q(p,:);
    % shape functions and first derivatives
    [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions    
    x        = N'*enode;
    plot(x(1),x(2),'black*');
  end
end
plot(nodes(13,1),nodes(13,2),'rs','MarkerSize',9,'MarkerFaceColor','red')
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
dtime  = 2*mesh.deltax/c;
dtime2 = dtime*dtime;
time   = 3;
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
 
  % Newton-Raphson iterations to solve for d(t+dtime)  
  disp ('=================================')
  disp (sprintf('   Time step %d',istep));
  disp ('=================================')
  disp ('  NR iter : L2-norm residual')
    
  error    = 1;
  iiter    = 0;
    
  while error > tol
    iiter    = iiter + 1;    
    nacce    = 1/(beta*dtime2)*(ndisp-dtilde);
    % compute the tangent and internal force
    [geoMat, stiffMat, fint] = computeTangentMatrix(mesh,material, ndisp);        
    % compute the external force
    [fext] = computeExternalForce(mesh,g,rho);
    % modified stiffness matrix
    K      = 1/(beta*dtime2)*massMat + stiffMat + geoMat;    
    % Residual vector
    res    = massMat*nacce + fint - fext;
    % Boundary conditions
    [K,res]  = applyDirichletBCs(K,res,udofs,vdofs,uFixed,vFixed);
    % Solving for displacement increment
    ddu      = -K\res;
    % Update displacements and do not forget to update the mesh
    ndisp     = ndisp + ddu;
    mesh.node = mesh.node + [ddu(1:2:dofCount) ddu(2:2:dofCount)];           
    error     = norm(ddu);
    disp (sprintf(' %s %i %s %5.4e ', 'Iter',iiter, ':', error) );
  end
 
  nacce    = 1/(beta*dtime2)*(ndisp-dtilde);
  nvelo    = nvelo0 + (1-gamma)*dtime*nacce0 + gamma*dtime*nacce;
  
  % advance to the next time step
  t        = t + dtime;
  istep    = istep + 1;
  
  nacce0   = nacce;
  nvelo0   = nvelo;
  ndisp0   = ndisp;
  
  % compute displacement at the GP of the marked element 
  esctr  = elements(markedElem,:);
  enode  = nodes(esctr,:);  
  edisp  = ndisp(2*esctr); % vertical displacements
  % second GP of the marked element
  p      = 2;
  pt   = Q(p,:);
  % shape functions and first derivatives
  [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
  uu       = N'*edisp;
  pDisp(istep) = uu;    % vertical displacement recorded 
  
  if (  mod(istep-1,interval) == 0 )   
    vtuFile = sprintf('%s%d','femImpVibratingCantilever',istep-1);
    VTKPostProcess(nodes,elements,2,'Quad4',vtuFile,sigma,...
      [ndisp(1:2:dofCount) ndisp(2:2:dofCount)]);
  end
end
%%
disp([num2str(toc),'   DONE ']);

ul = load('vibratingCantileverUL.mat');

figure
hold on
plot(ta,pDisp,'black-','LineWidth',1.8)
plot(ul.ta,ul.pDisp,'blue-','LineWidth',1.8)
xlabel('time')
ylabel('displacement')
set(gca,'FontSize',16)
legend('Implicit UL','Explicit UL')
grid on

%%
figure
fac=1e4;
plot_mesh(mesh.node,elements,elemType,'k-',1.);

