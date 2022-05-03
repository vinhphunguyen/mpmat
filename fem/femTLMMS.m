% This file implements the Total Lagrangian FEM.
% One dimensional problem with Method of Manufactured Solution (MMS).
% The grid consists of two-noded or three-noded linear elements.
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
addpath ('/Users/vingu/my-codes/matlab2tikz/src');

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
G      = 0.0001; %100*G
%G      = 0.05;  % large deformation

elemType  = 'L2';
%elemType  = 'L3';

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
%  Computational grid: linear or quadratic elements
L     = 1;
m     = [3 4 5 6 7 8 9 10 12]; % choose number of elements for convenrgence study
ne    = 2^m(9);

if strcmp(elemType,'L2')
 [mesh]=buildGrid1D(L,ne,0);
 noGPs = 2;
 mesh.deltax  = L/ne;
else
  mesh.nodeCount = 2*ne+1;
  mesh.node      = linspace(0,L,mesh.nodeCount)';
  mesh.elemCount = ne;
  elems          = zeros(ne,3);
  for e=1:ne
    elems(e,:) = 2*e-1:2*e+1;
  end
  elems(:,[2 3])   = elems(:,[3 2]);
  mesh.element = elems;
  mesh.deltax  = L/ne;
  noGPs = 3;
end

elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
elements  = mesh.element;
nodes     = mesh.node;             % to be updated in a UL formulation (Eulreian x)
deltax    = mesh.deltax;


%%
% Integration points:

[W,Q]   = quadrature( noGPs, 'GAUSS', 1 ); % 2 Gauss point rule
pCount  = noGPs*elemCount;             % total of Gauss points of the mesh
stress  = zeros(pCount,1);             % stress at all GPs
defGrad = ones(pCount,1);              % deformation gradient at all GPs
indices = zeros(elemCount,noGPs);

i = 1;
for e=1:elemCount
  for p=1:length(W)
    indices(e,p) = i; i = i + 1;
  end
end

%% Time loop
tol = 0;

dtime = 0.5*deltax/c;
time  = 0.02;
t     = 0.;
istep = 0;

nmass     = zeros(nodeCount,nodeCount);  % lumped mass matrix
nvelo     = zeros(nodeCount,1);  % nodal  velocity vector (final)
nvelo0    = zeros(nodeCount,1);  % nodal velocity vector (begin)
ndisp     = zeros(nodeCount,1);  % nodal displacement vector (end)
niforce   = zeros(nodeCount,1);  % nodal internal force vector
neforce   = zeros(nodeCount,1);  % nodal external force vector

% Compute lumped mass matrix once
% row-sum technique
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

% nodal integration technique

% gaussLobatto.W = [0.333333333333 0.333333333333 1.3333333333333];
% gaussLobatto.Q = [-1 1 0];
% 
% for e=1:elemCount
%   esctr = elements(e,:);
%   enode = nodes(esctr);
%   % loop over Gauss point
%   for p=1:length(gaussLobatto.W)
%     pt       = gaussLobatto.Q(p);
%     [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
%     J0       = enode'*dNdxi;                  % element Jacobian matrix
%     detJ     = det(J0);
%     mm       = N * N' * rho * detJ * gaussLobatto.W(p);
%     nmass(esctr,esctr) = nmass(esctr,esctr) + mm;
%   end
% end
% nmassd = diag(nmass);

% initialise nodal velocities
for i=1:nodeCount
  nvelo0(i) = mmsV(nodes(i),0);
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
    ue    = ndisp(esctr);  
    % loop over integration points
    for p=1:length(W)
      pt   = Q(p);
      % shape functions and first derivatives
      [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
      J0       = enode'*dNdxi;                  % element Jacobian matrix
      dNdx     = dNdxi/J0;
      wt       = W(p)*J0;
      xx       = dot(N,enode);                  % global Gauss point (for body force)
      pid      = indices(e,p);
      gradu    = dot(dNdx,ue);    % gradient velocity      
      F        = 1 + gradu;            
      defGrad(pid) = F;
      stress(pid)  = lambda*log(F)/F + mu*F - mu/F; % Neo-Hookean
      % internal force
      niforce(esctr) = niforce(esctr) - wt*stress(pid)*dNdx;      
      % external force due to manufactured body force
      neforce(esctr) = neforce(esctr) + rho*wt*N*mmsB(xx,t);      
    end
  end
  
  % update nodal velocity
  nforce    = niforce + neforce;
  acce      = nforce./nmassd;
  if (istep==0), acce = 0.5*acce; end
  nvelo     = nvelo0 + acce*dtime;
  nvelo(1)  = 0; nvelo(nodeCount)  = 0; % Boundary conditions f1 = m1*a1, a1=0
  deltaU    = dtime*nvelo;              % displacement increment
  ndisp     = ndisp + deltaU;           % displacement at the end of time step
  % update stresses at integration points
%   for e=1:elemCount
%     esctr = elements(e,:);
%     enode = nodes(esctr);
%     ue    = ndisp(esctr);   
%     % loop over integration points
%     for p=1:length(W)
%       pt   = Q(p);
%       % shape functions and first derivatives
%       [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
%       J0       = enode'*dNdxi;                  % element Jacobian matrix
%       dNdx     = dNdxi/J0;
%       wt       = W(p)*J0;
%       pid      = indices(e,p);
%       %v1       = nvelo(esctr(1));
%       %v2       = nvelo(esctr(2));      
%       %Lp      = dNdx(1) * v1 + dNdx(2) * v2;    % gradient velocity
%       %F       = (1 + Lp*dtime)*defGrad(pid);    
%       gradu    = dot(dNdx,ue);    % gradient velocity      
%       F        = 1+gradu;            
%       defGrad(pid) = F;
%       stress(pid)  = lambda*log(F)/F + mu*F - mu/F; % Neo-Hookean
%     end
%   end  
  % swap old and new velocities (this is not needed in MPM!!!)
  nvelo0 = nvelo;
    % advance to the next time step
  t = t + dtime;
  istep    = istep + 1;
  % compute displacement error norm
  dispNorm = 0;
  for e=1:elemCount
    esctr  = elements(e,:);
    enode  = nodes(esctr);
    edisp  = ndisp(esctr);
    % loop over integration points
    for p=1:length(W)
      pt   = Q(p);
      % shape functions and first derivatives
      [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
      J0       = enode'*dNdxi;                  % element Jacobian matrix
      %dNdx     = dNdxi/J0;
      wt       = W(p)*J0;
      X        = dot(N,enode);
      uh       = dot(N,edisp);
      uex      = mmsU(X,t);
      dispNorm = dispNorm + (uh-uex)^2*wt;
      %exact    = exact    + (uex)^2*wt;
    end
  end
 
  dispNorm = sqrt(dispNorm);
  %err      = max(err,dispNorm); % L_inf used in least-square PIC paper
  err(istep) = dispNorm;    % L2 norm used in CPDI paper  
  
  
end
%%
disp([num2str(toc),'   DONE ']);
%% convergence plot

% result for G=0.0001; small deformation
dispLinG1=[2.721087193526597e-06;
          6.839588526994928e-07;
          1.712134960211386e-07;
          4.281748314990817e-08;
          1.072139121871020e-08;
          2.678338662003873e-09;
          6.690890936948027e-10;
          1.673037509305871e-10;
          8.1683e-12
     ];
   
dispLinG2=[0.001355239633159;
           3.406751641746584e-04;
          8.528133903848277e-05;
          2.132730547758120e-05;
          5.339959396182416e-06;
          1.334058544737701e-06;
          3.332748019851188e-07;
          8.333394331780576e-08;
          5.2076e-09
     ];

size=[0.125000000000000;
      0.062500000000000;
      0.031250000000000;
      0.015625000000000;
      0.007812500000000;
      0.003906250000000;
      0.001953125000000;
      9.765625000000000e-04;
      2.4414e-04];


polyfit(log(size),log(dispLinG1),1)
polyfit(log(size),log(dispLinG2),1)

%plot(log(size),log(disp),'black*-')
% loglog(size,dispLinf1,'black*-','LineWidth',1.8)
% hold on
% loglog(size,dispLinf2,'reds-','LineWidth',1.8)
% loglog(size,dispLinf3,'cyans-','LineWidth',1.8)
% xlabel('Element side')
% ylabel('Error')
% set(gca,'FontSize',16)
% legend('small deformation','large deformation')
% grid on
%box on
%set(gca ,'YTickLabel',num2str(disp,1))
% axis([0 100 -0.15 0.2])

set(gca,'FontSize',16)

loglog(size,dispLinG1,'black*-','LineWidth',1.2)
hold on
loglog(size,dispLinG2,'reds-','LineWidth',1.2)
% loglog(size,dispLinf2,'cyans-','LineWidth',1.8)
xlabel('Element size')
ylabel('Error')
%set(gca,'FontSize',16)
legend('G=0.0001', 'G=0.05')
grid on
%%
figure
plot(ta(2:end),err);
set(gca,'FontSize',16)
%set(gca, 'YTickLabel', num2str(get(gca, 'YTick')))


