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
E      = 1000;               % Young modulus
rho    = 1000;              % density
c      = sqrt(E/rho);
v      = 0.1;

%% MLS weight functions
shape = 'circle' ;         % shape of domain of influence
dmax1 = 2.5 ;              % radius = dmax * nodal spacing
dmax2 = 2.5 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;   % using cubic spline weight function

%%
%  Computational grid: two-noded elements
L     = 1;
ne    = 20;
[mesh]=buildGrid1D(L,ne,0);

elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
elements  = mesh.element;
nodes     = mesh.node;
deltax    = mesh.deltax;
omegaC    = deltax;

%%
% Material points:
ppc       = 3; % # particles per cell
[mesh1]   = buildGrid1D(5*deltax,5,0);
particles = buildParticlesGeom(mesh1,ppc,rho);

pCount = 2*particles.pCount;
xp   = [particles.xp;particles.xp+L-5*deltax];
xp0  = xp;
vp   = [particles.vp;particles.vp];
Vp   = [particles.Vp;particles.Vp];
Vp0  = Vp;
Fp   = ones(pCount,1);
s    = [particles.s;particles.s];
eps  = [particles.eps;particles.eps];
Mp   = [particles.Mp;particles.Mp];
rhop = [particles.rho;particles.rho];

% initial velocities, initial stress=0
for p=1:pCount
  if xp(p) < L/2
    vp(p) = v;
  else
    vp(p) = -v;
  end
end
vp0=vp;

xp1 = xp(1:pCount/2);
xp2 = xp(pCount/2+1:end);

% Domain of influence for every nodes(sampling points) 
% Uniformly distributed nodes
% Definition : rad = dmax*deltaX

di1 = ones(length(xp),1)*dmax1*deltax;
di2 = ones(length(xp),1)*dmax2*deltax;
di3 = ones(length(xp1),1)*dmax2*deltax;

%%
figure
hold on
plot(nodes,zeros(nodeCount,1)+1/2,'black-s');
plot(xp1,zeros(pCount/2,1)+1/2,'b*');
plot(xp2,zeros(pCount/2,1)+1/2,'r*');
axis([0 L+0.1 0 1.])

%%
% data structure to store the material points for each element
% this data structure is updated for every time step

pElems  = ones(pCount ,1);
mpoints = cell (elemCount ,1);

for p=1:pCount
  x = xp(p);
  e = floor(x/deltax) + 1;
  pElems(p) = e; % particle "p" stays in element "e"
  for e=1:elemCount
    id = find(pElems==e);
    mpoints{e}=id ; % mpoints{e}?> indices of particles in "e"
  end
end

activeElems = unique(pElems);
activeNodes = unique(elements(activeElems,:));

%% nodal quantities
nmass       = zeros(nodeCount,1);  % nodal mass vector
nvelo       = zeros(nodeCount,1);  % nodal velocity vector (final)
nvelo0      = zeros(nodeCount,1);  % nodal velocity vector (begin)
niforce     = zeros(nodeCount,1);  % nodal internal force vector
neforce     = zeros(nodeCount,1);  % nodal external force vector

cellDensity = zeros(elemCount,1);  % cell-centerd density
cellStress  = zeros(elemCount,1);  % cell-centerd stress

for ie=1:length(activeElems)
    e     = activeElems(ie);
    esctr = elements(e,:);
    enode = nodes(esctr);        
    xc    = mean(enode);    
    index = defineSupport(xp,xc,di2);
    %if length(index) <= 2, disp('A singular'); end
    phi   = mlsLinearBasis1D(xc,index,xp,di2,form);
    cellDensity(e) = cellDensity(e) + dot(phi,rhop(index));    
 end


%% Time loop
tol = 0;

dtime = 0.001;
time  = 3.6;
t     = 0.;
istep = 0;

nsteps = floor(time/dtime);
err    = zeros(nsteps,1);

ta = []; ka = []; sa = []; pos={}; ix=1;
%%
while ( t < time )
  disp(['time step ',num2str(t)]);
  nmass(:)     = 0;
  nvelo0(:)    = 0;
  niforce(:)   = 0;
  neforce(:)   = 0;
  % loop over computational cells or elements
  for ie=1:length(activeElems)
    e     = activeElems(ie);
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
  end
  % project velocity to grid nodes (MLS)
  for ii=1:length(activeNodes)
    i     = activeNodes(ii);
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
  %nvelo(1)  = 0; 
  %nvelo(nodeCount)  = 0; % Boundary conditions
  
  % update particle velocity and position and stresses
  k = 0; u = 0;
  for ie=1:length(activeElems)
    e     = activeElems(ie);
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
      dEps    = dtime * Lp;
      s(pid)  = s(pid)   + E * dEps;
      eps(pid)= eps(pid) + dEps;
      %rhop(pid) = rho/F;
                  
      k = k + 0.5*vp(pid)^2*Mp(pid);
      u = u + 0.5*s(pid)*eps(pid)*Vp(pid);
    end
  end
  
  % update the element particle list
%   pe = floor(xp/deltax)+1;
%   for e=1:elemCount
%     id  = find(pe==e);
%     mpoints{e}=id;
%   end
  xp1 = xp(1:pCount/2);
  xp2 = xp(pCount/2+1:end);
  s1  = s(1:pCount/2);
  s2  = s(pCount/2+1:end);
  rho1  = rhop(1:pCount/2);
  rho2  = rhop(pCount/2+1:end);

  for p=1:pCount
    x = xp(p);
    e = floor(x/deltax) + 1;
    pElems(p) = e; % particle "p" stays in element "e"
    for e=1:elemCount
      id = find(pElems==e);
      mpoints{e}=id ; % mpoints{e}?> indices of particles in "e"
    end
  end
  
  activeElems = unique(pElems);
  activeNodes = unique(elements(activeElems,:));
    
  % project particle density/stress to grid centers (MLS)
  activeElems1 = unique(floor(xp1/deltax)+1);
  activeElems2 = unique(floor(xp2/deltax)+1);
  cellDensity(:) = 0;
  cellStress(:)  = 0;
  for ie=1:length(activeElems)
    e     = activeElems(ie);
    esctr = elements(e,:);
    enode = nodes(esctr);        
    xc    = mean(enode);    
    index = defineSupport(xp,xc,di2);
    %if length(index) <= 2, disp('A singular'); end
    phi   = mlsLinearBasis1D(xc,index,xp,di2,form);
    cellDensity(e) = cellDensity(e) + dot(phi,rhop(index));
    cellStress(e)  = cellStress(e)  + dot(phi,s(index));
  end
  
%   for ie=1:length(activeElems2)
%     e     = activeElems2(ie);
%     esctr = elements(e,:);
%     enode = nodes(esctr);        
%     xc    = mean(enode);    
%     index = defineSupport(xp2,xc,di3);
%     %if length(index) <= 2, disp('A singular'); end
%     phi   = mlsLinearBasis1D(xc,index,xp2,di3,form);
%     cellDensity(e) = cellDensity(e) + dot(phi,rho2(index));
%     cellStress(e)  = cellStress(e)  + dot(phi,s2(index));
%   end
  
  % advance to the next time step
  t = t + dtime;
  istep = istep + 1; 
  
  ta = [ta;t];
  ka = [ka;k];
  sa = [sa;u];
  
      if  ( mod(istep,100) == 0 )
    pos{ix} = xp; ix = ix+1;
    end
end
%%
disp([num2str(toc),'   DONE ']);


mpm=load('mpmCollision1D.mat');

figure 
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(mpm.ta(1:end),mpm.ka(1:end),'b--','LineWidth',1.6);
plot(ta(1:end),sa(1:end),'r-','LineWidth',2);
plot(mpm.ta(1:end),mpm.sa(1:end),'r--','LineWidth',2);
%plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy')
legend('kinetic-IMPM','kinetic-MPM','strain-IMPM','strain-MPM')
%%
% for i=1:length(pos)
%   figure
%   hold on
%   plot(nodes,zeros(nodeCount,1)+1/2,'r-s');
%   plot(pos{i},zeros(pCount,1)+1/2,'b*');
%   title( sprintf('%i',num2str(90*dtime*i) ));
%   frame(i) = getframe;
% end
% close all;
% movie(frame);
