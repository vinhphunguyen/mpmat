% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
%
% Shape functions: using standard FE shape function routine by converting
% particle position to natural coordinates.
%
% Spring mass system with standard GIMP.
%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 7 September 2015.

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
rho   = 0.1;         % density
mm    = 3.33;        % mass of the MASS
ms    = mm/10000;    % mass of the spring
A     = 0.1;
L     = 0.3;         % length of the spring
g     = 250;         % gravity acceleration
ks    = E*A/L;       % spring stiffness k
omega = sqrt(ks/mm); % natural frequency
Delta = mm*g/ks;     % spring stretch at static equilibrium


interval     = 100;% time interval for saving vtp files.
vtkFileName  = 'mpm2DTwoDisks';

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

l = A;
w = 0.4;

noX0      = 1;        % number of elements along X direction
noY0      = 4;        % number of elements along Y direction
ghostCell = 1;        % NOTE: for GIMP, ghostCell must be 1

[bGrid]   = buildGrid2D(l,w,noX0,noY0, ghostCell);

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

pCount = 10;
deltas = L/(pCount-1);
ratio  = deltay/deltas;
xx = 0 : deltas : L;

Mp  = zeros(pCount,1);                % mass
Fp  = repmat([1 0 0 1],pCount,1);     % gradient deformation
s   = zeros(pCount,3);                % stress
eps = zeros(pCount,3);                % strain
vp  = zeros(pCount,2);                % velocity
xp  = zeros(pCount,2);                % position
% pBasis = zeros(pCount,4);             % grid functions for particles
% pGradx = zeros(pCount,4);
% pGrady = zeros(pCount,4);


xp(:,1) = A/2 + deltax;
xp(:,2) = xx + w-L + deltay - deltas/2;

Mp(:)   = ms;
Mp(1)   = mm;
Vp      = Mp/rho;
Vp(1)   = mm/rho;
Vp0     = Vp;

y0      = xp(1,2);                    % store the initial pos. of the MASS

% GIMP particle size

lpx = deltax/2;
lpy = deltas;

%% find GIMP element connectivity

gimpElement = cell(elemCount,1);

for e=1:elemCount
    neighbors      = getNeighbors(e, bGrid.numx, bGrid.numy);
    neighborNodes  = element(neighbors,:);
    gimpElement{e} = unique(neighborNodes);
end

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

pElems  = ones(pCount,1);
mpoints = cell(elemCount,1);

for ip=1:pCount
  x = xp(ip,1); y = xp(ip,2);
  e = floor(x/bGrid.deltax) + 1 + bGrid.numx*floor(y/bGrid.deltay);
  pElems(ip) = e;
end

for ie=1:elemCount
  id  = find(pElems==ie);
  mpoints{ie}=id;
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
plot(xp(:,1),xp(:,2),'r.','markersize',20);
plot(xp(1,1),xp(1,2),'k.','markersize',40);
axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy
ya = [];           % displacement of the mass

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-12; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.4*bGrid.deltax/c;
time  = 0.1;        % simulation time
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < 0.05 )
  disp(['time step ',num2str(t)])
  % reset grid data
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
      % particle mass and momentum to node
      stress = s(pid,:);
      mp     = Mp(pid);
      px     = xp(pid,:);
      for i=1:length(esctr)
        id    = esctr(i);
        x     = px - node(id,:);
        [N,dNdx]=getGIMP2D(x,deltax,deltay,lpx,lpy);
        dNIdx = dNdx(1);
        dNIdy = dNdx(2);
        nmass(id)        = nmass(id)        + N*mp;
        nmomentum0(id,:) = nmomentum0(id,:) + N*mp*vp(pid,:);
        niforce(id,1)    = niforce(id,1) - Vp(pid)*(stress(1)*dNIdx + stress(3)*dNIdy);
        niforce(id,2)    = niforce(id,2) - Vp(pid)*(stress(3)*dNIdx + stress(2)*dNIdy);
        %if ( t == 0 )
          neforce(id,:)   = neforce(id,:) + N*mp*bodyF;
        %end
      end
    end
  end
  nforce                   = niforce + neforce;
  % update nodal momenta
  nmomentum                = nmomentum0 + nforce*dtime;
  nmomentum(bGrid.tNodes,:) = 0;
  
  % update particle velocity and position
  for e=1:elemCount
    esctr = element(e,:);    
    mpts  = mpoints{e};
    % loop over particles
    for p=1:length(mpts)
      pid  = mpts(p);      
      px   = xp(pid,:);   
      for i=1:length(esctr)
        id         = esctr(i);             
        x          = px - node(id,:);
        [N,~]      = getGIMP2D(x,deltax,deltay,lpx,lpy);
        massInv    = (1/nmass(id))*N;
        vp(pid,:)  = vp(pid,:) +  (nmomentum(id,:)-nmomentum0(id,:))*massInv;
        %xp(pid,:)  = xp(pid,:) + dtime * nmomentum(id,:)*massInv;                  
      end       
      % mapped back grid momenta (used to compute L,,epsilon and stress)
      for i=1:length(esctr)
        id         = esctr(i);                
        x          = px - node(id,:);
        [N,~]      = getGIMP2D(x,deltax,deltay,lpx,lpy);
        nmomentumS(id,:)  = nmomentumS(id,:) +  Mp(pid)*vp(pid,:)*N;                  
      end  
    end
  end
  %nmomentumS = nmomentum;
  nmomentumS(bGrid.tNodes,:) = 0;  % boundary conditions
  
  k = 0; u = 0;
  for e=1:elemCount
    esctr = element(e,:);    
    mpts  = mpoints{e};
    % loop over particles
    for p=1:length(mpts)
      pid  = mpts(p);
      px   = xp(pid,:);   
      Lp   = zeros(2,2);      
      for i=1:length(esctr)
        id = esctr(i);
        vI = [0 0];
        if nmass(id) > 0               
          x          = px - node(id,:);
          [N,dNdx]   = getGIMP2D(x,deltax,deltay,lpx,lpy);
          vI         = nmomentumS(id,:)/nmass(id);  % nodal velocity
          Lp         = Lp + vI'*dNdx;         % particle gradient velocity
          xp(pid,:)  = xp(pid,:) + dtime * vI * N;                  
        end        
      end
      
      F       = ([1 0;0 1] + Lp*dtime)*reshape(Fp(pid,:),2,2);
      Fp(pid,:)= reshape(F,1,4);
      Vp(pid) = det(F)*Vp0(pid);
      dEps    = dtime * 0.5 * (Lp+Lp');
      dsigma  = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
      s(pid,:)  = s(pid,:) + dsigma';
      eps(pid,:)= eps(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
      
      k = k + 0.5*(vp(pid,1)^2+vp(pid,2)^2)*Mp(pid);
      u = u + 0.5*Vp(pid)*s(pid,:)*eps(pid,:)';
    end
  end
  
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
  
  % store time,velocty for plotting
  ta = [ta;t];   ka = [ka;k]; sa = [sa;u]; ya = [ya;y0-xp(1,2)];
  
  % advance to the next time step
  t     = t + dtime;
  istep = istep + 1;
end

%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

t = 0:0.0005:time;
exact = (Delta)*(1-cos(omega*t));

figure
set(gca,'FontSize',14)
hold on
plot(t,exact,'bo','LineWidth',1.6);
plot(ta,ya,'r','LineWidth',1.6);
xlabel('Time')
ylabel('Displacement')
legend('exact','MPM')
%set(gca,'YTick',[0 0.001 0.002 0.003 0.004 0.005])
set(gca,'YTickLabel',num2str(get(gca,'YTick').'))
%axis([0 0.001 0 6])



% savefile = 'mpm-explicit-2disks-dt001.mat';
% save(savefile,'ta','ka','sa');

disp([num2str(toc),'   DONE '])
