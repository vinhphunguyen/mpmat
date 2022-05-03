% This file implements the Updated Lagrangian FEM.
% One dimensional problem with Method of Manufactured Solution (MMS).
% The grid is one two-noded linear element.
% Leapfrog time integration.
%
% This file contains commands to plot a convergence curve.
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 28 September 2015.

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
E      = 1e7;               % Young modulus
nu     = 0.0;               % Poisson ratio
rho    = 1000;              % density
lambda = E*nu/(1+nu)/(1-2*nu);
mu     = E/2/(1+nu);
c      = sqrt(E/rho);
G      = 0.0001; 
%G      = 100*G;
G      = 0.5;

% function handles for MMS (manufactured solutions)

mmsU = @(x,t)      G*sin(pi*x)*sin(c*pi*t);
mmsV = @(x,t) pi*c*G*sin(pi*x)*cos(c*pi*t);
mmsF = @(x,t) 1 + pi*G*cos(pi*x)*sin(c*pi*t);
mmsB = @(x,t) (1/rho)*pi^2*mmsU(x,t)*( (lambda/mmsF(x,t)/mmsF(x,t))*(1-log(mmsF(x,t))) + ...
  mu*(1+1/mmsF(x,t)/mmsF(x,t)) -E );

% [X,T] = meshgrid(0:0.05:1,0:0.0006:0.02);
% U     = G*sin(pi*X).*sin(c*pi*T);
% surf(X,T,U)
% xlabel('X')
% ylabel('t')
% zlabel('u')
% set(gca,'FontSize',16)
% %camlight right
% shading interp
% %axis equal

%%
%  Computational grid: two-noded elements
L     = 1;
m     = [3 4 5 6 7 8]; % choose number of elements for convenrgence study
ne    = 2^m(4);
[mesh]=buildGrid1D(L,ne,0);

elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
elements  = mesh.element;
nodes     = mesh.node;             % to be updated in a UL formulation (Eulreian x)
nodes0    = nodes;                 % material coordinates X
deltax    = mesh.deltax;
elemType  = 'L2';

%%
% Integration points:

[W,Q]   = quadrature( 2, 'GAUSS', 1 ); % 2 Gauss point rule
pCount  = length(W)*elemCount;         % total of Gauss points of the mesh
stress  = zeros(pCount,1);             % stress at all GPs
defGrad = ones(pCount,1);              % deformation gradient at all GPs
indices = zeros(elemCount,length(W));

i = 1;
for e=1:elemCount
  for p=1:length(W)
    indices(e,p) = i; i = i + 1;
  end
end

%% Time loop
tol = 0;

dtime = 0.2*deltax/c;
time  = 0.02;
t     = 0.;
istep = 0;

nmass     = zeros(nodeCount,nodeCount);  % lumped mass matrix
nvelo     = zeros(nodeCount,1);  % nodal  velocity vector (final)
ndisp     = zeros(nodeCount,1);  % nodal displacement vector (end)
niforce   = zeros(nodeCount,1);  % nodal internal force vector
neforce   = zeros(nodeCount,1);  % nodal external force vector

% Compute lumped mass matrix once
for e=1:elemCount
  esctr = elements(e,:);
  enode = nodes(esctr);
  % loop over Gauss point
  for p=1:length(W)
    pt       = Q(p);
    [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
    J0       = enode'*dNdxi;                  % element Jacobian matrix
    detJ     = det(J0);
    mm       = N * N' * rho * detJ * W(p);
    nmass(esctr,esctr) = nmass(esctr,esctr) + mm;
  end
end
% make the consistent mass matrix a diagonal (lumped) mass matrix
nmassd = sum(nmass,1)';

% initialise nodal velocities
for i=1:nodeCount
  nvelo(i)  = mmsV(nodes(i),0);
  ndisp(i)  = mmsU(nodes(i),0); % not neccessary in this case but
end

nsteps = floor(time/dtime);
err    = zeros(nsteps,1);
ta     = 0:dtime:time;

while ( t < time )
  disp(['time step ',num2str(t)]);
  niforce(:)   = 0;
  neforce(:)   = 0;
  % loop over computational cells or elements
  for e=1:elemCount
    esctr = elements(e,:);
    enode = nodes(esctr);
    enode0 = nodes0(esctr);
    ue    = ndisp(esctr);  
    % loop over integration points
    for p=1:length(W)
      pt   = Q(p);
      % shape functions and first derivatives
      [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
      J0       = enode'*dNdxi;                  % element Jacobian matrix
      dNdx     = dNdxi/J0;
      wt       = W(p)*J0;      
      pid      = indices(e,p);
      gradu    = dot(dNdx,ue);                  % grad of displacement
      F        = 1/(1-gradu);                  
      stress(pid)  = lambda*log(F)/F + mu*F - mu/F; % Neo-Hookean      
      % internal force
      niforce(esctr)    = niforce(esctr) - wt*stress(pid)*dNdx;      
      % external force due to manufactured body force
      XX              = dot(N,enode0);      % global Gauss point (for body force)
      neforce(esctr)  = neforce(esctr) + rho*W(p)*(enode0'*dNdxi)*N*mmsB(XX,t);      
    end
  end
  
  % update nodal velocity
  nforce    = niforce + neforce;
  acce      = nforce./nmassd;
  if (istep==0), acce = 0.5*acce; end
  nvelo     = nvelo + acce*dtime;
  nvelo(1)  = 0; nvelo(nodeCount)  = 0; % Boundary conditions f1 = m1*a1, a1=0
  deltaU    = dtime*nvelo;              % displacement increment
  ndisp     = ndisp + deltaU;           % displacement at the end of time step
    % advance to the next time step
  t = t + dtime;
  istep    = istep + 1;
  % update node coordinates (updated Lagrangian)  
  nodes = nodes + deltaU; 
  % compute displacement error norm
  dispNorm = 0;
  for e=1:elemCount
    esctr  = elements(e,:);
    enode  = nodes(esctr);
    enode0 = nodes0(esctr);
    edisp  = ndisp(esctr);
    % loop over integration points
    for p=1:length(W)
      pt   = Q(p);
      % shape functions and first derivatives
      [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
      J0       = enode'*dNdxi;                  % element Jacobian matrix
      %dNdx     = dNdxi/J0;
      wt       = W(p)*J0;
      X        = dot(N,enode0);
      uh       = dot(N,edisp);
      uex      = mmsU(X,t);
      dispNorm = dispNorm + (uh-uex)^2*wt;
      %exact    = exact    + (uex)^2*wt;
    end
  end 
  dispNorm = sqrt(dispNorm);  
  err(istep) = dispNorm;    % L2 norm used in CPDI paper      
end
%%
disp([num2str(toc),'   DONE ']);

%% convergence plot

% G=0.0001
dispLinf0=[2.721091084984657e-06;
          6.839626952106446e-07;
          1.712173043039249e-07;
          4.282127208716546e-08;
          1.072517079192418e-08;
          2.681292697981045e-09
     ];

% G=0.01;   
dispLinf1=[2.720662490135263e-04;
          6.838537185840609e-05;
          1.711870155361189e-05;
          4.281061852364640e-06;
          1.071939388213574e-06;
          2.679848470533935e-07
     ];   

% G=0.05;   
dispLinf2=[0.001355239633159;
          3.406751641746931e-04;
          8.528133903849269e-05;
          2.132730547761791e-05;
          5.339833481861222e-06;
          1.334042758337293e-06
     ];   
   
size=[0.125;
  0.062500000000000;
  0.031250000000000;
  0.015625000000000;
  0.007812500000000;
  0.003906250000000];

polyfit(log(size),log(dispLinf0),1)


loglog(size,dispLinf0,'reds-','LineWidth',1.8)
hold on
loglog(size,dispLinf1,'blues-','LineWidth',1.8)
loglog(size,dispLinf2,'blacks-','LineWidth',1.8)
xlabel('Element size')
ylabel('Error')
set(gca,'FontSize',16)
legend('G=0.0001','G=0.01', 'G=0.05')
grid on
%box on
%set(gca ,'YTickLabel',num2str(disp,1))
% axis([0 100 -0.15 0.2])
%%
figure
plot(ta(2:end),err);
set(gca,'FontSize',16)
%set(gca, 'YTickLabel', num2str(get(gca, 'YTick')))


