% This file implements the Improved Material Point Method.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
% Leapfrog time integration.
%
% Moving Least Square approximation used to construct particle data on the
% grid (nodes/centers).
% MLS introduces two parameters into the problem: smoothing length (domain of influence).
% Generally, vvelocity MLS has one smoothing length and density/stress MLS have
% another smoothing length. Parameter study is needed.
%
% Convergence test of MPM using the Method of Manufactured Solution (MMS).
% This file contains commands to plot a convergence curve.
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 9 October 2015.

%%

addpath ../grid/
addpath ../basis/
addpath ../particleGen/
addpath ../constitutiveModels/
addpath ../util/
addpath ../geoMesh/
addpath ../externals/
addpath ../postProcessing/
addpath ../mls/

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
G      = 0.1;

% function handles for MMS (manufactured solutions)

mmsU = @(x,t)      G*sin(pi*x)*sin(c*pi*t);
mmsV = @(x,t) pi*c*G*sin(pi*x)*cos(c*pi*t);
mmsF = @(x,t) 1 + pi*G*cos(pi*x)*sin(c*pi*t);
mmsB = @(x,t) (1/rho)*pi^2*mmsU(x,t)*( (lambda/mmsF(x,t)/mmsF(x,t))*(1-log(mmsF(x,t))) + ...
  mu*(1+1/mmsF(x,t)/mmsF(x,t)) -E );

%% MLS weight functions
shape = 'circle' ;         % shape of domain of influence
dmax1 = 2.0 ;              % radius = dmax * nodal spacing
dmax2 = 2.0 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;   % using cubic spline weight function

%%
%  Computational grid: two-noded elements
L     = 1;
m     = [3 4 5 6 7]; % choose number of elements for convenrgence study
ne    = 2^m(3);
[mesh]=buildGrid1D(L,ne,0);

elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
elements  = mesh.element;
nodes     = mesh.node;
deltax    = mesh.deltax;
omegaC    = deltax;

%%
% Material points:
ppc = 2; % # particles per cell
particles=buildParticlesGeom(mesh,ppc,rho);

pCount = particles.pCount;
xp   = particles.xp;
xp0  = particles.xp;
vp   = particles.vp;
Vp   = particles.Vp;
Vp0  = particles.Vp0;
Fp   = particles.Fp;
s    = particles.s;
eps  = particles.eps;
Mp   = particles.Mp;
rhop = particles.rho;
mbody = zeros(pCount,1);

% initial velocities, initial stress=0
for p=1:particles.pCount
  vp(p) = mmsV(xp(p),0);
end
vp0=vp;

% Domain of influence for every nodes(sampling points) 
% Uniformly distributed nodes
% Definition : rad = dmax*deltaX

di1 = ones(length(xp),1)*dmax1*deltax;
di2 = ones(length(xp),1)*dmax2*deltax;

%%
% hold on
% plot(nodes,zeros(nodeCount,1)+1/2,'r-s');
% plot(xp,zeros(particles.pCount,1)+1/2,'b*');
% axis([0 L+0.1 0 1.])

%%
% data structure to store the material points for each element
% this data structure is updated for every time step

pElems  = ones(particles.pCount ,1);
mpoints = cell (elemCount ,1);

for p=1:particles.pCount
  x = xp(p);
  e = floor(x/deltax) + 1;
  pElems(p) = e; % particle "p" stays in element "e"
  for e=1:elemCount
    id = find(pElems==e);
    mpoints{e}=id ; % mpoints{e}?> indices of particles in "e"
  end
  XX   = xp0(p);
  mbody(p) = mmsB(XX,0);
end

%% nodal quantities
nmass       = zeros(nodeCount,1);  % nodal mass vector
nvelo       = zeros(nodeCount,1);  % nodal velocity vector (final)
nvelo0      = zeros(nodeCount,1);  % nodal velocity vector (begin)
niforce     = zeros(nodeCount,1);  % nodal internal force vector
neforce     = zeros(nodeCount,1);  % nodal external force vector

cellDensity = zeros(elemCount,1);  % cell-centerd density
cellStress  = zeros(elemCount,1);  % cell-centerd stress
cellBody    = zeros(elemCount,1);  % cell-centerd body force
 
for e=1:elemCount
    esctr = elements(e,:);
    enode = nodes(esctr);        
    xc    = mean(enode);    
    index = defineSupport(xp,xc,di2);
    %if length(index) <= 2, disp('A singular'); end
    phi   = mlsLinearBasis1D(xc,index,xp,di2,form);
    cellDensity(e) = cellDensity(e) + dot(phi,rhop(index));
    cellStress(e)  = cellStress(e)  + dot(phi,s(index)); 
    cellBody(e)    = cellBody(e)    + dot(phi,mbody(index));
 end


%% Time loop
tol = 0;

dtime = 0.2*deltax/c;
time  = 0.02;
t     = 0.;
istep = 0;

nsteps = floor(time/dtime);
err    = zeros(nsteps,1);
ta     = 0:dtime:time;

while ( t < time )
  disp(['time step ',num2str(t)]);
  nmass(:)     = 0;
  nvelo0(:)    = 0;
  niforce(:)   = 0;
  neforce(:)   = 0;
  % loop over computational cells or elements
  for e=1:elemCount
    esctr = elements(e,:);
    enode = nodes(esctr);
    mpts  = mpoints{e};
    % one point quadrature
    xc   = mean(enode);    
    % shape functions and first derivatives
    N1  = 1 - abs(xc-enode(1))/deltax;
    N2  = 1 - abs(xc-enode(2))/deltax;
    dN1 = -1/deltax;
    dN2 =  1/deltax;
    % mass
    rhoC              = cellDensity(e);
    nmass(esctr(1))   = nmass(esctr(1))     + N1*rhoC*omegaC;
    nmass(esctr(2))   = nmass(esctr(2))     + N2*rhoC*omegaC;
    % internal force
    sigma             = cellStress(e);
    niforce(esctr(1)) = niforce(esctr(1)) - omegaC*sigma*dN1;
    niforce(esctr(2)) = niforce(esctr(2)) - omegaC*sigma*dN2;
    % external force due to manufactured body force    
    body              = cellBody(e);
    neforce(esctr(1)) = neforce(esctr(1)) + rhoC*omegaC*N1*body;
    neforce(esctr(2)) = neforce(esctr(2)) + rhoC*omegaC*N2*body;    
  end
  % project velocity to grid nodes (MLS)
  for i=1:nodeCount
    pt    = nodes(i);
    index = defineSupport(xp,pt,di1);
    %if length(index) <= 2, disp('A singular'); end
    phi   = mlsLinearBasis1D(pt,index,xp,di1,form);
    nvelo0(i) = nvelo0(i) + dot(phi,vp(index));
  end
  
  % update nodal velocity
  nforce        = niforce + neforce;
  % leapfrog integration scheme
  if ( istep == 0 ), nforce = 0.5*nforce; end
  nvelo     = nvelo0 + dtime*nforce./nmass;
  nvelo(1)  = 0; 
  nvelo(nodeCount)  = 0; % Boundary conditions
  % update particle velocity and position and stresses
  for e=1:elemCount
    esctr = elements(e,:);
    enode = nodes(esctr);
    mpts  = mpoints{e};
    % loop over particles
    for p=1:length(mpts)
      pid  = mpts(p);
      N1  = 1 - abs(xp(pid)-enode(1))/deltax;
      N2  = 1 - abs(xp(pid)-enode(2))/deltax;
      dN1 = -1/deltax;
      dN2 =  1/deltax;
      v1 = 0; v2 = 0;
      if nmass(esctr(1)) > tol
        vp(pid)  = vp(pid) + N1*(nvelo(esctr(1))-nvelo0(esctr(1)));
        v1       = nvelo(esctr(1));
        up       = dtime * N1*v1;
        xp(pid)  = xp(pid) + up;
      end
      
      if nmass(esctr(2)) > tol
        %vp(pid)  = vp(pid) + N2*nforce(esctr(2))/nmass(esctr(2));
        vp(pid)  = vp(pid) + N2*(nvelo(esctr(2))-nvelo0(esctr(2)));
        v2       = nvelo(esctr(2));
        up       = dtime * N2*v2;
        xp(pid)  = xp(pid) + up;
      end      
      % gradient velocity
      Lp      = dN1 * v1 + dN2 * v2;
      F       = (1 + Lp*dtime)*Fp(pid);
      Fp(pid) = F;
      Vp(pid) = F*Vp0(pid);
      s(pid)  = lambda*log(F)/F + mu*F - mu/F; % Neo-Hookean
      rhop(pid) = rho/F;
    end
  end
  
  % project particle density/stress to grid centers (MLS)
  cellDensity(:) = 0;
  cellStress(:)  = 0;
  for e=1:elemCount
    esctr = elements(e,:);
    enode = nodes(esctr);        
    xc    = mean(enode);    
    index = defineSupport(xp,xc,di2);
    %if length(index) <= 2, disp('A singular'); end
    phi   = mlsLinearBasis1D(xc,index,xp,di2,form);
    cellDensity(e) = cellDensity(e) + dot(phi,rhop(index));
    cellStress(e)  = cellStress(e)  + dot(phi,s(index));
  end
  
  % update the element particle list
  pe = floor(xp/deltax)+1;
  for e=1:elemCount
    id  = find(pe==e);
    mpoints{e}=id;
  end
  
  % advance to the next time step
  t = t + dtime;
  istep = istep + 1;
  
  % compute error norm
  dispNorm = 0;
  for p=1:pCount
    xx0 = xp0(p);
    xx  = xp(p);
    up  = xx - xx0;
    uex = mmsU(xx0,t);
    dispNorm = dispNorm + Vp(p)*(up-uex)^2;
    
    XX   = xp0(p);
    mbody(p) = mmsB(XX,t);  
  end
  dispNorm = sqrt(dispNorm);
  err(istep) = dispNorm;
  
end
%%
disp([num2str(toc),'   DONE ']);


%% convergence plot

% G=0.0001
disp1=[1.181767235896483e-05;
       3.062114801402464e-06;
       7.712732858447164e-07;
       1.931993470892981e-07];

% G=0.1
disp2=[0.013097808721005;
       0.007842436073303;
       0.008800404624961;
  2.758304861721332e-04];     
 
% G=0.001, standard MPM     
disp1MPM=[1.078575308827472e-05;
      2.563045809675697e-06;
      6.089991540886196e-07;
      1.478400791855386e-07];     


size=[0.125000000000000;
  0.062500000000000;
  0.031250000000000;
  0.015625000000000];

polyfit(log(size),log(disp1),1)
polyfit(log(size),log(disp2),1)

loglog(size,disp1,'black*-','LineWidth',1.8)
hold on
loglog(size,disp1MPM,'blacks-','LineWidth',1.8)
loglog(size,disp2,'red*-','LineWidth',1.8)
xlabel('Element size')
ylabel('Error')
%legend('MPM','Exact')
set(gca,'FontSize',16)
grid on
legend('G=0.0001-IMPM','G=0.0001-MPM','G=0.01-IMPM','G=0.01-MPM')
%set(gca ,'YTickLabel',num2str(disp,1))
% axis([0 100 -0.15 0.2])

%%
figure
plot(ta(2:end),err);
set(gca,'FontSize',16)


