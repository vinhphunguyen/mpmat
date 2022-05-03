% This file implements the Generalized Interpolator Material Point Method.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Particles are generated using an unstructured FE mesh (gmsh).
%
%
% Two elastic disks come into contact.
% Example taken from Sulsky et al. 1994 paper.
%
% Fast implementation using C-MEX functions.
% NOTE:
% 1. bodies = findActiveElemsAndNodes(bodies,mesh) works only for MPM.
% 2. GIMP: nmass(activeNodes) sometimes zero at extra nodes!!!
%
% 14 August:
%  two disks stick together.
%  If comment in MEX file defo[ip]=fxx etc. or comment vol[ip]=vol0[ip*detF 
%  then it works.
%
% Vinh Phu Nguyen
% University of Adelaide, Adelaide, Australia.
% August 2014.

%%

addpath ../fem_util/
addpath ../fem-functions/
addpath ../post-processing/
addpath mex/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%

E   = 1000;        % Young's modulus
nu  = 0.3;         % Poisson ratio
rho = 1000;        % density
kappa = 3-4*nu;    % Kolosov constant
mu    = E/2/(1+nu);% shear modulus
v     = 0.1;   % initial particle velocity

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell= 1;
lx       = 1;
ly       = 1;
numx2    = 30;       % number of elements along X direction
numy2    = 30;       % number of elements along Y direction
[mesh]   = buildGrid2D(lx,ly,numx2,numy2, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;

%% particle generation
% body1

ppc    = 2; % # of particle per cell is ppc x ppc
[W,Q]=quadrature(  ppc, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

dx = mesh.deltax/(ppc);
dy = mesh.deltay/(ppc);

if (ghostCell)
  mesh.lpx = dx;
  mesh.lpy = dy;
end

center = [0.2 0.2];
radius = 0.2;

volume = []; mass   = []; coord  = [];

for e=1:elemCount                 % start of element loop
  sctr = element(e,:);          %  element scatter vector
  pts  = node(sctr,:);
  x1 = pts(1,:); % first corner of the cell
  for i=1:ppc
    for j=1:ppc
      x(1) = x1(1) + dx*0.5 + (j-1)*dx;
      x(2) = x1(2) + dy*0.5 + (i-1)*dy;
      r  = norm(x-center);
      if ( r-radius < 0 )
        volume  = [volume;dx*dy];
        mass    = [mass; dx*dy*rho];
        coord   = [coord;x];
      end
    end
  end
end

% for e=1:elemCount                 % start of element loop
%     sctr = element(e,:);          %  element scatter vector
%     pts  = node(sctr,:);
%
%     for q=1:size(W,1)                           % quadrature loop
%         pt=Q(q,:);                              % quadrature point
%         wt=W(q);                                % quadrature weight
%         [N,dNdxi]=lagrange_basis('Q4',pt);
%         J0 = pts'*dNdxi;
%         x  = N'*pts;
%         r  = norm(x-center);
%         if ( r-radius < 0 )
%             volume  = [volume;wt*det(J0)];
%             mass    = [mass; wt*det(J0)*rho];
%             coord   = [coord;x];
%         end
%     end
% end

coord(:,1) = coord(:,1) + mesh.deltax;
coord(:,2) = coord(:,2) + mesh.deltay;

pCount = length(volume);

% stored in body1 structure

body1.coord   = coord;
body1.volume  = volume;
body1.volume0 = volume;
body1.mass    = mass;
body1.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress  = zeros(pCount,3);                % stress
body1.strain  = zeros(pCount,3);                % strain
body1.velo    = ones(pCount,2)*v;               % velocity
body1.C       = C;
%body1.deform0 = body1.deform;
body1.gravity = 0;

% body 2
center = [0.8 0.8];
radius = 0.2;

volume = []; mass   = []; coord  = [];

for e=1:elemCount                 % start of element loop
  sctr = element(e,:);            %  element scatter vector
  pts  = node(sctr,:);
  
  x1 = pts(1,:); % first corner of the cell
  for i=1:ppc
    for j=1:ppc
      x(1) = x1(1) + dx*0.5 + (j-1)*dx;
      x(2) = x1(2) + dy*0.5 + (i-1)*dy;
      r  = norm(x-center);
      if ( r-radius < 0 )
        volume  = [volume;dx*dy];
        mass    = [mass; dx*dy*rho];
        coord   = [coord;x];
      end
    end
  end
end

% for e=1:elemCount                 % start of element loop
%     sctr = element(e,:);          %  element scatter vector
%     pts  = node(sctr,:);
%
%     for q=1:size(W,1)                           % quadrature loop
%         pt=Q(q,:);                              % quadrature point
%         wt=W(q);                                % quadrature weight
%         [N,dNdxi]=lagrange_basis('Q4',pt);
%         J0 = pts'*dNdxi;
%         x  = N'*pts;
%         r  = norm(x-center);
%         if ( r-radius < 0 )
%             volume  = [volume;wt*det(J0)];
%             mass    = [mass; wt*det(J0)*rho];
%             coord   = [coord;x];
%         end
%     end
% end

coord(:,1) = coord(:,1) + mesh.deltax;
coord(:,2) = coord(:,2) + mesh.deltay;

body2.volume  = volume;
body2.volume0 = volume;
body2.mass    = mass;
body2.coord   = coord;
body2.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body2.stress  = zeros(pCount,3);                % stress
body2.strain  = zeros(pCount,3);                % strain
body2.velo    = -ones(pCount,2)*v;              % velocity
body2.C       = C;
%body2.deform0 = body2.deform;
body2.gravity = 0;

bodies    = cell(2,1);
bodies{1} = body1;
bodies{2} = body2;
bodyCount = length(bodies);

%% find GIMP element connectivity

gimpElement = cell(elemCount,1);

for e=1:elemCount
  neighbors      = getNeighbors(e, mesh.numx, mesh.numy);
  neighborNodes  = element(neighbors,:);
  gimpElement{e} = unique(neighborNodes);
end

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

% The folowing function is only for MPM
%bodies = findActiveElemsAndNodes(bodies,mesh);

for ib=1:length(bodies)
  body      = bodies{ib};
  coord     = body.coord;
  elems     = ones(size(coord,1),1);
  
  for ip=1:size(coord,1)
    x = coord(ip,1); y = coord(ip,2);
    e = floor(x/mesh.deltax) + 1 + mesh.numx*floor(y/mesh.deltay);
    elems(ip) = e;
  end
  
  bodies{ib}.elements = unique(elems);
  bodies{ib}.nodes    = unique(vertcat(gimpElement{bodies{ib}.elements}));
  
  mpoints = cell(elemCount,1);
  for ie=1:elemCount
    id  = find(elems==ie);
    mpoints{ie}=id;
  end
  
  bodies{ib}.mpoints  = mpoints;
end


%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,2);  % nodal momentum vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector (no need for this exam)
nvelo     = zeros(nodeCount,2);  % nodal velocity vector
nacce     = zeros(nodeCount,2);  % nodal acceleration vector

%% plot mesh, particles

figure
set(gca,'FontSize',14)
hold on
plot_mesh(node,element,'Q4','k-',1.6);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','cy-',1.6);
xp1 = bodies{1}.coord;
xp2 = bodies{2}.coord;
plot(xp1(:,1),xp1(:,2),'k.','markersize',15);
plot(xp2(:,1),xp2(:,2),'r.','markersize',15);
plot(node(bodies{1}.nodes,1),node(bodies{1}.nodes,2),'*');
plot(node(bodies{2}.nodes,1),node(bodies{2}.nodes,2),'*');

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.001;
%time  = 60; % low velocity, used to test implicit
time  = 3.5;% explicit, high velocity
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
  disp(['time step ',num2str(t)])
  
  %reset grid data
%     nmass(:)     = 0;
%     nmomentum(:) = 0;
%     niforce(:)   = 0;
%     for ib=1:bodyCount                %% loop over bodies
%       body      = bodies{ib};
%       elems     = body.elements;
%       mpoints   = body.mpoints;
%       for ie=1:length(elems)         % loop over computational cells or elements
%         e     = elems(ie);
%         esctr = gimpElement{e};    % GIMP (extended) element connectivity        
%         mpts  = mpoints{e};        % particles inside element e
%         for p=1:length(mpts)       % loop over particles
%           pid    = mpts(p);
%           xp     = body.coord(pid,:);
%           stress = body.stress(pid,:);
%           Mp     = body.mass(pid);
%           vp     = body.velo(pid,:);
%           Vp     = body.volume(pid);
%           for i=1:length(esctr)  % loop over nodes of particle p
%             id    = esctr(i);
%             x     = xp - node(id,:);
%             [N,dNdx]=getGIMP2D(x,mesh.deltax,mesh.deltay,dx,dy);
%             dNIdx = dNdx(1);
%             dNIdy = dNdx(2);
%             nmass(id)       = nmass(id)       + N*Mp;
%             nmomentum(id,:) = nmomentum(id,:) + N*Mp*vp;
%             niforce(id,1)   = niforce(id,1) - Vp*(stress(1)*dNIdx + stress(3)*dNIdy);
%             niforce(id,2)   = niforce(id,2) - Vp*(stress(3)*dNIdx + stress(2)*dNIdy);
%           end
%         end
%       end
%     end
  
  % MEX function
  [nmass,nmomentum,niforce] = ParticlesToNodesGIMP(bodies,mesh);
  
  % update nodal momenta
  
  activeNodes=[bodies{1}.nodes; bodies{2}.nodes];
  nmomentum(activeNodes,:) = nmomentum(activeNodes,:) + niforce(activeNodes,:)*dtime;
  
  % compute updated velocities and accelerations
  nvelo(activeNodes,1) = nmomentum(activeNodes,1)./nmass(activeNodes);
  nvelo(activeNodes,2) = nmomentum(activeNodes,2)./nmass(activeNodes);
  
  nacce(activeNodes,1) = niforce(activeNodes,1)./nmass(activeNodes);
  nacce(activeNodes,2) = niforce(activeNodes,2)./nmass(activeNodes);
  
  % be careful with extra nodes in GIMP, where nmass=0 as N=0.
  idx=find(isnan(nvelo(:,1)));  
  nvelo(idx,:) = zeros(length(idx),2);
  nacce(idx,:) = zeros(length(idx),2);
  
  % note that there is no boundary conditions for this example
  
  UpdateParticlesGIMP(bodies,mesh,nvelo,nacce,dtime); % MEX function
  
%     for ib=1:bodyCount
%       body      = bodies{ib};
%       elems     = body.elements;
%       mpoints   = body.mpoints;
%       % loop over computational cells or elements
%       for ie=1:length(elems)
%         e     = elems(ie);
%         esctr = gimpElement{e};    % GIMP (extended) element connectivity
%         mpts  = mpoints{e};        % particles inside element e
%         for p=1:length(mpts)       % loop over particles
%           pid  = mpts(p);
%           xp   = body.coord(pid,:);
%           Mp   = body.mass(pid);
%           vp   = body.velo(pid,:);
%           Vp   = body.volume(pid);
%           Lp   = zeros(2,2);
%           for i=1:length(esctr)
%             id = esctr(i);
%             vI = nvelo(id,:);
%             aI = nacce(id,:);
%             x     = xp - node(id,:);
%             [N,dNdx]=getGIMP2D(x,mesh.deltax,mesh.deltay,dx,dy);
%             vp  = vp  + dtime * N*aI;
%             xp  = xp  + dtime * N*vI;
%             Lp  = Lp + vI'*dNdx;
%           end
%   
%           bodies{ib}.velo(pid,:) = vp;
%           bodies{ib}.coord(pid,:)= xp;
%   
%           % update stress last
%           F       = ([1 0;0 1] + Lp*dtime)*reshape(bodies{ib}.deform(pid,:),2,2);
%           bodies{ib}.deform(pid,:) = reshape(F,1,4);
%           bodies{ib}.volume(pid  ) = det(F)*bodies{ib}.volume0(pid);
%           dEps    = dtime * 0.5 * (Lp+Lp');
%           dStrain = dtime*[Lp(1) Lp(4) Lp(2)+Lp(3)]; % Voigt notation
%           [stress,strain] = updateStressHooke(dStrain,C,bodies{ib}.stress(pid,:),bodies{ib}.strain(pid,:));
%   
%           bodies{ib}.strain(pid,:)  = strain;
%           bodies{ib}.stress(pid,:)  = stress;
%   
%         end
%       end
%     end
  
  %   bodies{1}.coord0 = bodies{1}.coord;
  %   bodies{2}.coord0 = bodies{2}.coord;
  
  % compute kinetic and strain energies
  k = 0; u = 0;
  for ib=1:bodyCount
    body      = bodies{ib};
    for p=1:length(body.mass)
      vp   = body.velo(p,:);
      k = k + 0.5*(vp(1)^2+vp(2)^2)*body.mass(p);
      u = u + 0.5*body.volume(p)*body.stress(p,:)*body.strain(p,:)';
    end
    %bodies{ib}.deform0=bodies{ib}.deform;
  end
  
  % update the element particle list
  
  for ib=1:length(bodies)
    body      = bodies{ib};
    coord     = body.coord;
    elems     = ones(size(coord,1),1);
    
    for ip=1:size(coord,1)
      x = coord(ip,1); y = coord(ip,2);
      e = floor(x/mesh.deltax) + 1 + mesh.numx*floor(y/mesh.deltay);
      elems(ip) = e;
    end
    
    bodies{ib}.elements = unique(elems);
    bodies{ib}.nodes    = unique(vertcat(gimpElement{bodies{ib}.elements}));
    
    mpoints = cell(elemCount,1);
    for ie=1:elemCount
      id  = find(elems==ie);
      mpoints{ie}=id;
    end
    
    bodies{ib}.mpoints  = mpoints;
  end
  
  % store time,velocty for plotting
  
  ta = [ta;t];
  ka = [ka;k];
  sa = [sa;u];
  
  % advance to the next time step
  
  t     = t + dtime;
  istep = istep + 1;
end

%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

%ka = ka*1000;
%sa = sa*1000;

ss=load('gimp-2disks-energies.mat');

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(ta(1:end),sa(1:end),'r-','LineWidth',2);
%plot(ss.ta(1:end),ss.ka(1:end),'b.','LineWidth',1.6);
%plot(ss.ta(1:end),ss.sa(1:end),'r.','LineWidth',2);
plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy (x1E-3)')
legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%axis([0 3 0 3])

disp([num2str(toc),'   DONE '])
