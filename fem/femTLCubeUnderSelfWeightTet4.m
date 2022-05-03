% This file implements the Total Lagrangian FEM.
% The grid consists of four-noded tetrahedral elements.
% Leapfrog time integration.
%
% Large deformation of a cube under self weight.
%
% Vinh Phu Nguyen
% Monash University
% 30 May 2019.

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
nu     = 0.3;               % Poisson ratio
rho    = 1050e-12;              % density
lambda = E*nu/(1+nu)/(1-2*nu);
mu     = E/2/(1+nu);
c      = sqrt(E/rho);
g      = -1e6;

identity = [1 0 0;0 1 0; 0 0 1]; % identity matrix
nsd      = 3; % number of spatial dimension

%% Computational mesh
meshFile = 'cube.msh';
mesh     = load_gmsh (meshFile);

elemType   = 'H4';
nodeCount  = mesh.nbNod;
elemCount  = mesh.nbTets;
nodes      = mesh.POS(:,1:3);
elements   = mesh.TETS(1:elemCount,1:4);
nodes0     = nodes;

%markedNode1 = find(abs(nodes(:,1)-4)<1e-10);
%markedNode2 = find(abs(nodes(:,2)-b)<1e-10);
%markedNode3 = find(abs(nodes(:,3)-0)<1e-10);

markedNode  = 1;%intersect(intersect(markedNode1,markedNode2),markedNode3);

fixNodes=find(abs(nodes(:,2)-1000)<1e-10);

%%
% Integration points:
noGPs   = 1;
[W,Q]   = quadrature( 1, 'TRIANGULAR', 3 ); % 1 Gauss point rule
pCount  = noGPs*elemCount;             % total of Gauss points of the mesh
stress  = zeros(pCount,2,2);           % stress at all GPs, 2x2 matrix for 1 GP
indices = zeros(elemCount,noGPs);

i = 1;
for e=1:elemCount
  for p=1:length(W)
    indices(e,p) = i; i = i + 1;
  end
end

hemin = 1e6;

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
    J0       = enode'*dNdxi;                  % element Jacobian matrix
    detJ     = det(J0);
    vol      = W(p)*detJ;
    he       = nthroot(vol,3);
    hemin    = min(he,hemin);
  end
end
plot(nodes(13,1),nodes(13,2),'rs','MarkerSize',9,'MarkerFaceColor','red')
axis off


%% Time loop
dtime = 0.2*hemin/c;
time  = 0.25;
t     = 0.;
istep = 0;
interval = 10;

nmass     = zeros(nodeCount,nodeCount);  % lumped mass matrix
nacce     = zeros(nodeCount,nsd);  % nodal  acceleration
nvelo     = zeros(nodeCount,nsd);  % nodal  velocity vector 
ndisp     = zeros(nodeCount,nsd);  % nodal displacement vector
niforce   = zeros(nodeCount,nsd);  % nodal internal force vector
neforce   = zeros(nodeCount,nsd);  % nodal external force vector

sigma     = zeros(nodeCount,3);  % nodal external force vector

%Compute lumped mass matrix once
%row-sum technique
for e=1:elemCount
  esctr = elements(e,:);
  enode = nodes(esctr,:);
  % loop over Gauss point
  for p=1:length(W)
    pt       = Q(p,:);
    [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
    J0       = enode'*dNdxi;                  % element Jacobian matrix
    detJ     = det(J0);
    mm       = N * N' * rho * detJ * W(p);
    nmass(esctr,esctr) = nmass(esctr,esctr) + mm;
  end
end
% make the consistent mass matrix a diagonal (lumped) mass matrix
nmassd = 1./sum(nmass,1)'; % already inversed

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


while ( t < time )
  disp(['time step ',num2str(t)]);
  niforce(:)   = 0;
  neforce(:)   = 0;
  % loop over elements
  for e=1:elemCount
    esctr  = elements(e,:);
    enode  = nodes(esctr,:);    
    enode0 = nodes0(esctr,:);    
    ue     = ndisp(esctr,:)';  
    % loop over integration points
    for p=1:length(W)
      pt       = Q(p,:);
      % shape functions and first derivatives
      [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
      J0       = enode'*dNdxi;                  % element Jacobian matrix
      dNdx     = dNdxi/J0;                      % equals B0^T 
      wt       = W(p)*det(J0);      
      pid      = indices(e,p);
      F        = identity + ue*dNdx; 
      %F        = enode'*dNdx; 
      invF     = inv(F);
      detF     = det(F);
      %if detF < 0, error('negative J'); end
      P        = mu*invF*(F*F'-identity) + lambda*log(detF)*invF;
      %stress(pid,:,:)  = P;
      % internal force
      niforce(esctr,:) = niforce(esctr,:) - wt*dNdx*P;            
      % external force due to gravity
      neforce(esctr,2) = neforce(esctr,2) + rho*wt*N*g;                  
    end
  end
  
  % update nodal velocity
  nforce    = niforce + neforce;
  acce(:,1) = nforce(:,1).*nmassd;
  acce(:,2) = nforce(:,2).*nmassd;
  acce(:,3) = nforce(:,3).*nmassd;
  if (istep==0), acce = 0.5*acce; end
  nvelo     = nvelo + acce*dtime;
  % boundary conditions, left edge fixed
  nvelo(fixNodes,1)  = 0.;  
  nvelo(fixNodes,2)  = 0.; 
  % update displacement and mesh
  deltaU    = dtime*nvelo;              % displacement increment
  ndisp     = ndisp + deltaU;           % displacement at the end of time step 
  % advance to the next time step
  t        = t + dtime;
  istep    = istep + 1;
  
  % compute displacement at the marked node
  pDisp(istep) = ndisp(markedNode,2);    % vertical displacement recorded 
  
%   if (  mod(istep-1,interval) == 0 )   
%     vtuFile = sprintf('%s%d','femULVibratingCantilever',istep-1);
%     VTKPostProcess(nodes,elements,2,'Quad4',vtuFile,sigma,ndisp);
%   end
end
%%
disp([num2str(toc),'   DONE ']);

ss = load('../example3D/cpdi/cpdiTet4Cantilever.mat');

figure
hold on
plot(ss.ta,ss.ka,'red-','LineWidth',1.8)
plot(ta(1:end),pDisp,'blue-','LineWidth',1.8)
xlabel('time')
ylabel('displacement')
set(gca,'FontSize',16)
legend('CPDI-Tet4','FEM')
grid on

save('femTet4Cantilever1000.mat', 'ta', 'pDisp');

%%
figure
fac=1;
plot_mesh(nodes+fac*ndisp,elements,elemType,'k-',1.);
