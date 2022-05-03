% This file implements the Total Lagrangian FEM.
% The grid consists of four-noded quadrilateral elements.
% Leapfrog time integration.
%
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 2 October 2015.

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
nu     = 0.3;               % Poisson ratio
rho    = 1000;              % density
lambda = E*nu/(1+nu)/(1-2*nu);
mu     = E/2/(1+nu);
c      = sqrt(E/rho);
G      = 0.0001; %100*G
G      = 500*G;

elemType  = 'Q4';
identity  = [1 0;0 1]; 

% Manufactured solutions
pi1    = 0;%(2./3.)*pi;
mmsU1  = @(x,t)      G*sin(pi*x)*sin(c*pi*t);
mmsU2  = @(y,t)      G*sin(pi*y)*sin(c*pi*t+pi);

mmsV1  = @(x,t)  pi*c*G*sin(pi*x)*cos(c*pi*t);
mmsV2  = @(y,t)  pi*c*G*sin(pi*y)*cos(c*pi*t+pi);

mmsF11 = @(x,t) 1 + pi*G*cos(pi*x)*sin(c*pi*t);
mmsF22 = @(y,t) 1 + pi*G*cos(pi*y)*sin(c*pi*t+pi);

mmsB1  = @(x,y,t) (1/rho)*pi^2*mmsU1(x,t)*( (lambda*(1-log( mmsF11(x,t)*mmsF22(y,t) ))+mu)/mmsF11(x,t)/mmsF11(x,t) + ...
  mu - E );
mmsB2  = @(x,y,t) (1/rho)*pi^2*mmsU2(y,t)*( (lambda*(1-log( mmsF11(x,t)*mmsF22(y,t) ))+mu)/mmsF22(y,t)/mmsF22(y,t) + ...
  mu - E );

%% Computational grid
l = 1.;
w = 1.;

m         = [3 4 5 6 7]; % choose number of elements for convenrgence study
noX0      = 2^m(4);      % number of elements along X direction
noY0      = noX0;        % number of elements along Y direction
ghostCell = 0;
[bGrid]   = buildGrid2D(l,w,noX0,noY0, ghostCell);
nodes     = bGrid.node;
elements  = bGrid.element;
elemCount = bGrid.elemCount;
nodeCount = bGrid.nodeCount;

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

%% Time loop
dtime = 0.2*bGrid.deltax/c;
time  = 0.02;
t     = 0.;
istep = 0;
interval = 10;

nmass     = zeros(nodeCount,nodeCount);  % lumped mass matrix
nacce     = zeros(nodeCount,2);  % nodal  acceleration
nvelo     = zeros(nodeCount,2);  % nodal  velocity vector 
ndisp     = zeros(nodeCount,2);  % nodal displacement vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector

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
for i=1:nodeCount
  nvelo(i,1) = mmsV1(nodes(i,1),0);
  nvelo(i,2) = mmsV2(nodes(i,2),0);
  
  ndisp(i,1) = mmsU1(nodes(i,1),0);
  ndisp(i,2) = mmsU2(nodes(i,2),0);
end

nsteps = floor(time/dtime);
err    = zeros(nsteps,1);
ta     = 0:dtime:time;

while ( t < time )
  disp(['time step ',num2str(t)]);
  niforce(:)   = 0;
  neforce(:)   = 0;
  % loop over elements
  for e=1:elemCount
    esctr = elements(e,:);
    enode = nodes(esctr,:);    
    ue    = ndisp(esctr,:)';  
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
      % external force due to manufactured body force
      xx  = N'*enode;
      neforce(esctr,1) = neforce(esctr,1) + rho*wt*N*mmsB1(xx(1),xx(2),t);      
      neforce(esctr,2) = neforce(esctr,2) + rho*wt*N*mmsB2(xx(1),xx(2),t);      
    end
  end
  
  % update nodal velocity
  nforce    = niforce + neforce;
  acce(:,1) = nforce(:,1).*nmassd;
  acce(:,2) = nforce(:,2).*nmassd;
  if (istep==0), acce = 0.5*acce; end
  nvelo     = nvelo + acce*dtime;
  % boundary conditions
  nvelo(bGrid.lNodes,1)  = 0.;  
  nvelo(bGrid.rNodes,1)  = 0.; 
  nvelo(bGrid.bNodes,2)  = 0.;  
  nvelo(bGrid.tNodes,2)  = 0.;  
  deltaU    = dtime*nvelo;              % displacement increment
  ndisp     = ndisp + deltaU;           % displacement at the end of time step  
  
  % advance to the next time step
  t        = t + dtime;
  istep    = istep + 1;
  
  % compute displacement error norm
  dispNorm = 0; 
  for e=1:elemCount
    esctr = elements(e,:);
    enode = nodes(esctr,:);    
    edisp = ndisp(esctr,:);  
    % loop over integration points
    for p=1:length(W)
      pt   = Q(p,:);
      % shape functions and first derivatives
      [N,dNdxi]= lagrange_basis(elemType,pt);   % element shape functions
      J0       = enode'*dNdxi;                  % element Jacobian matrix      
      wt       = W(p)*det(J0);      
      xx       = N'*enode;     
      uu       = N'*edisp; 
      uex1     = mmsU1(xx(1),t);
      uex2     = mmsU2(xx(2),t);
      dispNorm = dispNorm + ((uu(1)-uex1)^2 + (uu(2)-uex2)^2)*wt;
    end
  end
  dispNorm = sqrt(dispNorm);
  %err      = max(err,dispNorm); % L_inf used in least-square PIC paper
  err(istep) = dispNorm;    % L2 norm used in CPDI paper  
  
%   if (  mod(istep-1,interval) == 0 )   
%     vtuFile = sprintf('%s%d','mmsAxisAligned',istep-1);
%     VTKPostProcess(nodes,elements,2,'Quad4',vtuFile,sigma,ndisp);
%   end
end
%%
disp([num2str(toc),'   DONE ']);

%% convergence plot

% result for G=0.0001; small deformation
dispLinf1=[2.394868713559189e-06;
           5.921706068451673e-07;
           1.476413273705449e-07;
          3.747990103299733e-08          
     ];
   
% result for G=0.05, large deformation, same time step   
dispLinf2=[0.001268873324708;
           3.161080534381289e-04;
           7.895706013747854e-05;
          1.973399433197621e-05
     ];   
   


hh=[0.125;
  0.062500000000000;
  0.031250000000000;
  0.015625000000000];


polyfit(log(hh),log(dispLinf1),1)


loglog(hh,dispLinf1,'reds-','LineWidth',1.8)
hold on
loglog(hh,dispLinf2,'cyans-','LineWidth',1.8)
xlabel('Element size')
ylabel('Error')
set(gca,'FontSize',16)
legend('small deformation','large deformation')
grid on
%%
figure
plot(ta(2:end),err);
set(gca,'FontSize',16)
%set(gca, 'YTickLabel', num2str(get(gca, 'YTick')))

%%
figure
fac=6;
plot_mesh(nodes+fac*ndisp,elements,elemType,'k-',1.);
%%

tt  = 0:0.001:0.02;
fxx = G*sin(pi*(1/6))*sin(c*pi*tt);
fyy = G*sin(pi*(1/3))*sin(c*pi*tt+pi);

figure
hold on
plot(tt,fxx,'red','LineWidth',2)
plot(tt,fyy,'blue','LineWidth',2)
legend('ux','uy')
set(gca,'FontSize',16)
xlabel('time')
ylabel('displacement')
grid on

%plot(tt,fxx.*fyy,'cyan')
