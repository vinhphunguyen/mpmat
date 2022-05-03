% This file implements the Updated Lagrangian FEM.
% The grid consists of four-noded quadrilateral elements.
% Leapfrog time integration.
%
% Generalised vortex ring problem.
%
% Vinh Phu Nguyen
% Monash, Australia
% 6 January 2015.

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
E   = 1000;    % Young's modulus
nu  = 0.3;     % Poisson ratio
rho = 1000;    % density
K   = E/3/(1-2*nu);    % bulk modulus
mu    = E/2/(1+nu);% shear modulus
lambda = K - 2/3*mu;


I  = [1 0;0 1];

%% Computational grid
meshFile = 'ring.msh';
pmesh     = load_gmsh (meshFile);

elemType   = 'Q4';
nodeCount  = pmesh.nbNod;
elemCount  = pmesh.nbQuads;
nodes      = pmesh.POS(:,1:2);
elements   = pmesh.QUADS(1:elemCount,1:4);
nodes0     = nodes;

%%
% Integration points:
noGPs   = 4;
[W,Q]   = quadrature( 2, 'GAUSS', 2 ); % 2x2 Gauss point rule
pCount  = noGPs*elemCount;             % total of Gauss points of the mesh
stress  = zeros(pCount,2,2);           % stress at all GPs, 2x2 matrix for 1 GP
indices = zeros(elemCount,noGPs);

i = 1;
for e=1:elemCount
  for p=1:length(W)
    indices(e,p) = i; i = i + 1;
  end
end

markedElem = 6;

figure('Color','white')
hold on
plot_mesh(nodes,elements,elemType,'k-',1.);
plot_mesh(nodes,elements(markedElem,:),elemType,'r-',1.);
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
axis off


%% Time loop
c     = sqrt(E/rho);
time  = 1;
t     = 0.;
istep = 0;
interval = 2;                    % output write
intervalRemesh = 10;             % remeshing interval

nmass     = zeros(nodeCount,nodeCount);  % lumped mass matrix
nacce     = zeros(nodeCount,2);  % nodal  acceleration
nvelo     = zeros(nodeCount,2);  % nodal  velocity vector 
ndisp     = zeros(nodeCount,2);  % nodal displacement vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector

sigma     = zeros(nodeCount,3);  % nodal external force vector

%Compute lumped mass matrix once
%row-sum technique
h = 1e10;
for e=1:elemCount
  esctr = elements(e,:);
  enode = nodes(esctr,:);
  % loop over Gauss point
  area = 0;
  for p=1:length(W)
    pt       = Q(p,:);
    [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
    J0       = enode'*dNdxi;                  % element Jacobian matrix
    detJ     = det(J0);
    mm       = N * N' * rho * detJ * W(p);
    nmass(esctr,esctr) = nmass(esctr,esctr) + mm;
    area = area + detJ * W(p);
  end
  h = min(h,area);
end
% make the consistent mass matrix a diagonal (lumped) mass matrix
nmassd = 1./sum(nmass,1)'; % already inversed

dtime = 0.1*sqrt(h)/c;

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
    %enode0 = nodes0(esctr,:);    
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
      F        = inv(I - ue*dNdx); 
      %F        = enode'*dNdx; 
      invF     = inv(F);
      detF     = det(F);
      %if detF < 0, error('negative J'); end
      sig      = (mu*(F*F'-I) + lambda*log(detF)*I)/detF;
      %stress(pid,:,:)  = P;
      % internal force
      niforce(esctr,:) = niforce(esctr,:) - wt*dNdx*sig;      
      % external force due to MMS body force
      xg               = N'*enode;
      [bx,by]          = vortexBodyForces(xg(1),xg(2),t,mu,rho,G,Ri,Ro);
      neforce(esctr,:) = neforce(esctr,:) + rho*wt*N*[bx by];                  
    end
  end
  
  % update nodal velocity
  nforce    = niforce + neforce;
  acce(:,1) = nforce(:,1).*nmassd;
  acce(:,2) = nforce(:,2).*nmassd;
  if (istep==0), acce = 0.5*acce; end
  nvelo     = nvelo + acce*dtime;
  % boundary conditions, left edge fixed  
  %nvelo(bGrid.tNodes,2)  = 0.; 
  % update displacement and mesh
  deltaU    = dtime*nvelo;              % displacement increment
  ndisp     = ndisp + deltaU;           % displacement at the end of time step 
  nodes     = nodes + deltaU;           % updated node positions
  % advance to the next time step
  t        = t + dtime;
  istep    = istep + 1;
  
  if (  mod(istep,interval) == 0 )   
    vtuFile = sprintf('%s%d','femULVortex',istep-1);
    VTKPostProcess(nodes,elements,2,'Quad4',vtuFile,sigma,ndisp);
  end
end
%%
disp([num2str(toc),'   DONE ']);

%%
figure
fac=1;
plot_mesh(nodes,elements,elemType,'k-',1.);
axis on
%%


