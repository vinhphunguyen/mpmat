% This file implements the Total Lagrangian FEM.
% The grid consists of four-noded quadrilateral elements.
% Leapfrog time integration.
%
%
% Vinh Phu Nguyen
% 2019

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
E   = 1e6;    % Young's modulus
nu  = 0.3;     % Poisson ratio
rho = 1050;    % density
K   = E/3/(1-2*nu);    % bulk modulus
mu    = E/2/(1+nu);% shear modulus
lambda = K - 2/3*mu;
c      = sqrt(E/rho);

G      = 1;
Ri = 0.75;
Ro = 1.25;
T  = 1;
Rb = 0.5*(Ri+Ro);
Rm = Ro-Ri;

identity  = [1 0;0 1];

% folder for VTU files
folder = 'vtug02';
[s,d,a]= mkdir(folder);

%% Computational grid
meshFile = 'ring.msh';
pmesh     = load_gmsh (meshFile);

elemType   = 'Q4';
nodeCount  = pmesh.nbNod;
elemCount  = pmesh.nbQuads;
nodes      = pmesh.POS(:,1:2);
elements   = pmesh.QUADS(1:elemCount,1:4);
nodes0     = nodes;

innerNodes = [];
outerNodes = [];

epsilon = 1e-10;
for i = 1: length(nodes)
    x = nodes(i,1);
    y = nodes(i,2);
    radius  = sqrt ( x*x + y*y );
    if      ( abs (radius - Ri) < epsilon )
        innerNodes = [innerNodes;i];
    end    
    if ( abs (radius - Ro) < epsilon )
        outerNodes = [outerNodes;i];
    end
end

figure('Color','white')
hold on
plot_mesh(nodes0,elements,elemType,'k-',1.);
plot(nodes(innerNodes,1),nodes(innerNodes,2),'rs','MarkerSize',9,'MarkerFaceColor','red')
plot(nodes(outerNodes,1),nodes(outerNodes,2),'rs','MarkerSize',9,'MarkerFaceColor','cyan')
axis off


%%
% Integration points:
noGPs   = 4;
[W,Q]   = quadrature( 2, 'GAUSS', 2 ); % 2x2 Gauss point rule
pCount  = noGPs*elemCount;             % total of Gauss points of the mesh
stress  = zeros(pCount,2,2);           % stress at all GPs, 2x2 matrix for 1 GP
indices = zeros(elemCount,noGPs);

%% Time loop

time  = 1;
t     = 0.;
istep = 0;
interval = 100;

nmass     = zeros(nodeCount,nodeCount);  % lumped mass matrix
nacce     = zeros(nodeCount,2);  % nodal  acceleration
nvelo     = zeros(nodeCount,2);  % nodal  velocity vector 
ndisp     = zeros(nodeCount,2);  % nodal displacement vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector

sigma     = zeros(nodeCount,3);  % nodal external force vector

hemin = 1e12;
%Compute lumped mass matrix once
%row-sum technique
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
  hemin = min(hemin,area);
end
% make the consistent mass matrix a diagonal (lumped) mass matrix
nmassd = 1./sum(nmass,1)'; % already inversed


fac   = 0.8;
dtime = fac*hemin/c;

nsteps = floor(time/dtime);
err    = zeros(nsteps,1);
ta     = 0:dtime:time;

while ( t < time )
  disp(['time step ',num2str([istep t])]);
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
      if detF < 0, error('negative J'); end
      %P        = mu*invF*(F*F'-identity) + lambda*log(detF)*invF;
      P        = mu*(F-invF') + lambda*log(detF)*invF';
      %stress(pid,:,:)  = P;
      % internal force
      niforce(esctr,:) = niforce(esctr,:) - wt*dNdx*P';      
      % external force due to manufactured body force
%       xx  = N'*enode;
%       neforce(esctr,1) = neforce(esctr,1) + rho*wt*N*mmsB1(xx(1),xx(2),t);      
%       neforce(esctr,2) = neforce(esctr,2) + rho*wt*N*mmsB2(xx(1),xx(2),t);  
      xg  = N'*enode;
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
  
  % boundary conditions
  nvelo(innerNodes,1)  = 0.;  
  nvelo(innerNodes,2)  = 0.; 
  nvelo(outerNodes,1)  = 0.;  
  nvelo(outerNodes,2)  = 0.;  
  
  deltaU    = dtime*nvelo;              % displacement increment
  ndisp     = ndisp + deltaU;           % displacement at the end of time step  
  nodes0    = nodes0 + ndisp;

% exact solution

%   gt = G*sin(pi*t/T);
%   
%   for i = 1: length(nodes)
%     x  = nodes(i,1);
%     y  = nodes(i,2);
%     theta0 = atan2(y,x)/pi;
%     sigma(i,1) = theta0;
%     sigma(i,2) = theta0;
%     sigma(i,3) = theta0;
%     x0 = [x;y];
%     radius  = sqrt ( x*x + y*y );
%     hR  = 1 - 8 * (radius-Rb)^2/Rm^2 + 16*(radius-Rb)^4/Rm^4;
%     alpha = gt*hR;
%     Q     = [cos(alpha) -sin(alpha);...
%              sin(alpha) cos(alpha)];
%     x     = Q*x0;
%     u     = x - x0;
%     ndisp(i,1) = u(1);
%     ndisp(i,2) = u(2);
%   end

  
  
  % advance to the next time step
  t        = t + dtime;
  istep    = istep + 1;
  

  if (  mod(istep-1,interval) == 0 )   
    vtuFile = sprintf('%s%d','vtug02/femTLVortexG02',istep-1);
    VTKPostProcess(nodes0,elements,2,'Quad4',vtuFile,sigma,ndisp);
  end
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
fac=0.1;
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
