% This file implements the Material Point Method of Sulsky 1994 coupled
% with FEM for 1D membrane problem.
%
% The bGrid is a structured mesh consisting of 4-noded bilinear elements (Q4).
%
%
% Shape functions: using standard FE shape function routine by converting
% particle position to natural coordinates.
%
% Spring mass system with coupled FEM-MPM.
%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 10 September 2015.

%%

addpath ../../grid/
addpath ../../basis/
addpath ../../particleGen/
addpath ../../constitutiveModels/
addpath ../../util/
addpath ../../geoMesh/
addpath ../../externals/PolyMesher/
addpath ../../externals/
addpath ../../postProcessing/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%

E     = 1e6;         % Young's modulus
nu    = 0.0;         % Poisson ratio
rho   = 1000;        % density of MASS
mm    = 3.33;        % mass of the MASS
ms    = mm/10000;    % mass of the spring
A     = 0.1;
L     = 0.3;         % length of the spring
g     = 250;         % gravity acceleration

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);

option = 0; % option={0,1} to differentiate two methods of computing element forces

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational bGrid

l = A;
w = 0.4;

noX0      = 1;        % number of elements along X direction
noY0      = 8;        % number of elements along Y direction
ghostCell = 0;

[bGrid]    = buildGrid2D(l,w,noX0,noY0, ghostCell);

node      = bGrid.node;
element   = bGrid.element;
deltax    = bGrid.deltax;
deltay    = bGrid.deltay;
elemCount = bGrid.elemCount;
nodeCount = bGrid.nodeCount;
numx2     = bGrid.numx;
numy2     = bGrid.numy;
Vc        = deltax * deltay;

%% particles

pCount     = 20;
[pGrid]    = buildGrid1D(L,pCount-1, ghostCell);

Mp  = zeros(pCount,1);                % mass
Vp  = zeros(pCount,1);                % mass
Fp  = repmat([1 0 0 1],pCount,1);     % gradient deformation
s   = zeros(pCount,3);                % stress
eps = zeros(pCount,3);                % strain
vp  = zeros(pCount,2);                % velocity
xp  = zeros(pCount,2);                % position
fp  = zeros(pCount,2);                % particle forces (due to FEM-MPM)
up  = zeros(pCount,2);                % particle displacement (due to FEM-MPM)
le  = zeros(pCount-1,1);              % element length (due to FEM-MPM)
stress = zeros(pCount-1,1);           % element stress (due to FEM-MPM)
strain = zeros(pCount-1,1);           % element stress (due to FEM-MPM)
pBasis = zeros(pCount,4);             % bGrid functions for particles
pGradx = zeros(pCount,4);
pGrady = zeros(pCount,4);

le(:)   = pGrid.deltax;
xp(:,1) = A/2;
xp(:,2) = pGrid.node + w-L - 0.001;

Mp(:)   = ms;
Mp(1)   = mm;
Mp(end) = ms;
Vp(:)   = A*pGrid.deltax;
Vp0     = Vp;

y0      = xp(1,2);                    % store the initial pos. of the MASS
xp0 = xp;

L0 = L;%+ pGrid.deltax;
ks    = E*A/L0;       % spring stiffness k
omega = sqrt(ks/mm);  % natural frequency
Delta = mm*g/ks;      % spring stretch at static equilibrium
rhos  = ms/A/L0;      % density of spring (ideal: massless spring)


%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

pElems  = ones(pCount,1);
mpoints = cell(elemCount,1);

for p=1:pCount
  x = xp(p,1);
  y = xp(p,2);
  e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
  pElems(p) = e;
end

for e=1:elemCount
  id  = find(pElems==e);
  mpoints{e}=id;
end

%% node quantities

nmass      = zeros(nodeCount,1);  % nodal mass vector
nmomentum  = zeros(nodeCount,2);  % nodal momentum vector
nmomentum0 = zeros(nodeCount,2);  % nodal momentum vector
nmomentumS = zeros(nodeCount,2);  % nodal momentum vector (mapped back)
niforce    = zeros(nodeCount,2);  % nodal internal force vector
neforce    = zeros(nodeCount,2);  % nodal external force vector (no need for this exam)
bodyF      = [0 -g];

%% plot mesh, particles

figure(1)
hold on
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(xp(:,1),xp(:,2),'k.','markersize',20);
plot(xp(1,1),xp(1,2),'k.','markersize',40);
axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy
ya = [];           % displacement of the mass

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-12; % mass tolerance

c      = sqrt(E/rhos);
dtime1 = 0.8*bGrid.deltax/c;
dtime2 = 0.8*pGrid.deltax/c;
dtime  = min(dtime1,dtime1);
time   = 0.02;        % simulation time
t      = 0;
k      = 0; u = 0; ym = 0;

nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;
%%
while ( t < time )
  disp(['time step ',num2str(t)])
  % reset bGrid data
  nmass(:)      = 0;
  nmomentum0(:) = 0;
  nmomentumS(:) = 0;
  niforce(:)    = 0;
  neforce(:)    = 0;
  
  % loop over computational cells or elements
  for e=1:elemCount
    esctr = element(e,:);      % element connectivity
    enode = node(esctr,:);     % element node coords
    mpts  = mpoints{e};        % particles inside element e
    % loop over particles
    for p=1:length(mpts)
      pid  = mpts(p);
      pt(1)= (2*xp(pid,1)-(enode(1,1)+enode(2,1)))/deltax;
      pt(2)= (2*xp(pid,2)-(enode(2,2)+enode(3,2)))/deltay;
      
      [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
      J0       = enode'*dNdxi;             % element Jacobian matrix
      invJ0    = inv(J0);
      dNdx     = dNdxi*invJ0;
      pBasis(pid,:) = N;
      pGradx(pid,:) = dNdx(:,1);
      pGrady(pid,:) = dNdx(:,2);
      % particle mass and momentum to node
      
      mp     = Mp(pid);
      for i=1:length(esctr)
        id    = esctr(i);
        dNIdx = dNdx(i,1);
        dNIdy = dNdx(i,2);
        nmass(id)        = nmass(id)        + N(i)*mp;
        nmomentum0(id,:) = nmomentum0(id,:) + N(i)*mp*vp(pid,:);
        niforce(id,:)    = niforce(id,:)    + N(i)*fp(pid,:);
        neforce(id,:)    = neforce(id,:)    + N(i)*mp*bodyF;        
      end
    end
  end
  nforce                   = niforce + neforce;
  % update nodal momenta
  nmomentum                = nmomentum0 + nforce*dtime;
  nmomentum(bGrid.tNodes,:) = 0;
  %nmomentum0(bGrid.tNodes,:) = 0;
  
  % store time,velocty for plotting
  ta = [ta;t];   ka = [ka;k]; sa = [sa;u]; ya = [ya;ym];
  
  % update particle velocity and position
  for e=1:elemCount
    esctr = element(e,:);    
    mpts  = mpoints{e};
    % loop over particles
    for p=1:length(mpts)
      pid  = mpts(p);      
      N    = pBasis(pid,:);      
      for i=1:length(esctr)
        id         = esctr(i);                
        massInv    = (1/nmass(id))*N(i);
        vp(pid,:)  = vp(pid,:) +  (nmomentum(id,:)-nmomentum0(id,:))*massInv;
        xp(pid,:)  = xp(pid,:) + dtime * nmomentum(id,:)*massInv;                  
      end       
      % mapped back bGrid momenta (used to compute L,,epsilon and stress)
      for i=1:length(esctr)
        id = esctr(i);                                
        nmomentumS(id,:)  = nmomentumS(id,:) +  Mp(pid)*vp(pid,:)*N(i);                  
      end  
    end
  end

  if sum(isnan(vp(:,2))) 
    error('NAN')
  end
  
  nmomentumS(bGrid.tNodes,:) = 0;  % boundary conditions
  up(:) = 0;
  for e=1:elemCount
    esctr = element(e,:);    
    mpts  = mpoints{e};
    % loop over particles
    for p=1:length(mpts)
      pid  = mpts(p);
      Lp   = zeros(2,2);      
      N    = pBasis(pid,:);      
      for i=1:length(esctr)
        id = esctr(i);
        vI = [0 0];
        if nmass(id) > 0                    
          vI  = nmomentumS(id,:)/nmass(id);  % nodal velocity
        end
        up(pid,:) = up(pid,:) + dtime*vI*N(i); % particle increment displacement
      end
    end
  end
  
  fp(:) = 0; u = 0;
  
  if (option==1)
    for e=1:size(pGrid.element,1)
      esctr    = pGrid.element(e,:);
      nodes    = xp(esctr,:);
      ll       = norm(nodes(1,:)-nodes(2,:));
      cosTheta = (nodes(2,1)-nodes(1,1))/ll;
      sinTheta = (nodes(2,2)-nodes(1,2))/ll;
      Q        = [cosTheta  sinTheta;-sinTheta cosTheta];
      Qinv     = [cosTheta -sinTheta;sinTheta cosTheta];
      incDisp  = up(esctr,:);
      deltaU1  = Q*incDisp(1,:)';
      deltaU2  = Q*incDisp(2,:)';            
      incStrain = (1/ll)*(-deltaU1(1)+deltaU2(1));      
      %if incStrain < 0, incStrain=0.; end
      %if incStrain < 0, disp('negative inc. strain'); end
      strain(e) = strain(e) + incStrain;
      sig       = stress(e) + E*incStrain;
      %if strain(e) < 0, disp('NEGATIVE STRAIN');sig=0;end      
      %if sig < 0, sig = 0; end
      stress(e) = sig;
      le(e)     = ll;
      force     = sig*A;

      fp(esctr(1),:)  = fp(esctr(1),:) + (Qinv*[force;0])';
      fp(esctr(2),:)  = fp(esctr(2),:) - (Qinv*[force;0])';
      
      u = u + 0.5*le(e)*stress(e)*strain(e)*A;
    end
  else        
    % strain computed based on length, (also working well)
    for e=1:size(pGrid.element,1)
      esctr = pGrid.element(e,:);
      nodes = xp(esctr,:);
      ll    = norm(nodes(1,:)-nodes(2,:));
      epsil = (ll - le(e))/le(e);
      strain(e) = strain(e) + epsil;
      %if epsil < 0, epsil=0.; end
      sig   = stress(e) + E*epsil;
      stress(e) = sig;
      le(e)     = ll;
      force     = sig*A;
      cosTheta = (nodes(2,1)-nodes(1,1))/ll;
      sinTheta = (nodes(2,2)-nodes(1,2))/ll;
      fp(esctr(1),:)  = fp(esctr(1),:) + force*[ cosTheta sinTheta];
      fp(esctr(2),:)  = fp(esctr(2),:) - force*[ cosTheta sinTheta];
      
      u = u + 0.5*ll*sig*strain(e)*A;
    end
  end
  ym = y0-xp(1,2);            
  k  = 0.5*dot((vp(:,2).^2),Mp);
  
  % store time,velocty for plotting
  pos{istep} = xp;
  vel{istep} = vp;
  
  % update the element particle list
  for p=1:pCount
    x = xp(p,1);
    y = xp(p,2);
    e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
    pElems(p) = e;
  end
  
  for e=1:elemCount
    id  = find(pElems==e);
    mpoints{e}=id;
  end
  
  % advance to the next time step
  t     = t + dtime;
  istep = istep + 1;
  xp0 = xp;
end

%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

% mm=1.665;
% omega = sqrt(ks/mm); % natural frequency
% Delta = mm*g/ks;     % spring stretch at static equilibrium


tt     = 0:0.0001:time;
exact  = (Delta)*(1-cos(omega*tt));               % displaecment
ke     = 0.5*mm*(omega*Delta)^2*sin(omega*tt).^2; % kinetic energy
pe     = 0.5*ks*exact.^2;                         % potential energy

figure
set(gca,'FontSize',14)
hold on
plot(tt,exact,'bo','LineWidth',1.6);
plot(ta,ya,'r','LineWidth',1.6);
xlabel('Time')
ylabel('Displacement')
legend('exact','MPM')
%set(gca,'YTick',[0 0.001 0.002 0.003 0.004 0.005])
set(gca,'YTickLabel',num2str(get(gca,'YTick').'))
grid on
box on
%axis([0 0.001 0 6])
%% plot energies

figure
set(gca,'FontSize',14)
hold on
plot(tt,ke,'bo','LineWidth',1.6);
plot(tt,pe,'rs','LineWidth',1.6);
plot(ta,ka,'b-','LineWidth',1.6);
plot(ta,sa,'r-','LineWidth',1.6);
xlabel('Time')
ylabel('Energy')
legend('KE-exact','PE-exact','KE-MPM','PE-MPM')
%set(gca,'YTick',[0 0.001 0.002 0.003 0.004 0.005])
set(gca,'YTickLabel',num2str(get(gca,'YTick').'))
grid on, box on
%%
figure
hold on
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(xp0(:,1),xp0(:,2),'ro','markersize',20);
plot(xp(:,1),xp(:,2),'k.','markersize',20);
axis off

%%
ss = load('spring-mass-mpm.mat');

% figure
% set(gca,'FontSize',14)
% hold on
% plot(ss.ta,ss.ya,'b','LineWidth',1.6);
% plot(ta,ya,'r','LineWidth',1.6);
% xlabel('Time')
% ylabel('Displacement')
% legend('MPM','FEM-MPM')
% %set(gca,'YTick',[0 0.001 0.002 0.003 0.004 0.005])
% set(gca,'YTickLabel',num2str(get(gca,'YTick').'))
% grid on
% box on

% savefile = 'mpm-explicit-2disks-dt001.mat';
% save(savefile,'ta','ka','sa');

disp([num2str(toc),'   DONE '])
