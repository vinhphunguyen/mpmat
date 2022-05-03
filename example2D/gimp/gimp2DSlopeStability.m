% This file implements the Generalized Material Point Method (GIMP).
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Particles are generated as either Gauss points of the background mesh or
% as points that regularly divide the grid cells.
%
% Slope stability
% Soil: Drucker-Prager enhanced with ...
%
%
% Vinh Phu Nguyen
% 13 August 2014, Adelaide, Australia.

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

useMex = 0;

%% Material properties
%

E   = 20;               % Young's modulus
nu  = 0.42;             % Poisson ratio
rho = 1e-9;             % density
c   = 1e-3;             % cohesion
phi = 42*pi/180;        % angle of friction
psi = 0;                % angle of dilatation
mu  = 0.6;              % coefficient of friction (Coulomb friction)
g   = 9.81e3;           % gravitation
K0  = nu/(1-nu);
v0  = 1;                % imposed velocity

vtkFileName  = 'gimp2DSlopeStability';
interval     = 50;

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell= 1;
lx       = 22e3;
ly       = 11e3;
numx2    = 44;       % number of elements along X direction
numy2    = 22;      % number of elements along Y direction
[mesh]   = buildGrid2D(lx,ly,numx2,numy2, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;

%% particle generation
% slope

xA = [2e3  0];
xD = [12e3 10e3];

ppc    = 3; % # of particle per cell is ppc x ppc
useGPs = 1; % particles at Gauss points
useGPs = 0; % particles at geometric points regularly divide the elements

if (useGPs==0)
  dx = mesh.deltax/(ppc);
  dy = mesh.deltay/(ppc);
end

[W,Q]=quadrature(  3, 'GAUSS', 2 ); % 2x2 Gaussian quadrature


volume = []; mass   = []; coord  = [];

for e=1:elemCount                 % start of element loop
  sctr = element(e,:);          %  element scatter vector
  pts  = node(sctr,:);
  if (useGPs)
    for q=1:size(W,1)                           % quadrature loop
      pt=Q(q,:);                              % quadrature point
      wt=W(q);                                % quadrature weight
      [N,dNdxi]=lagrange_basis('Q4',pt);
      J0 = pts'*dNdxi;
      x  = N'*pts;
      if ( x(1) <= 4e3 )
        volume  = [volume;wt*det(J0)];
        mass    = [mass; wt*det(J0)*rho];
        coord   = [coord;x];
      end
    end
  else
    x1 = pts(1,:); % first corner of the cell
    for i=1:ppc
      for j=1:ppc
        x(1) = x1(1) + dx*0.5 + (j-1)*dx;
        x(2) = x1(2) + dy*0.5 + (i-1)*dy;
        if ( x(2) <= 10e3 )
          delta   = (xA(1)-x(1))*(xD(2)-x(2))-(xD(1)-x(1))*(xA(2)-x(2));
          if (delta < -1e-6)
            volume  = [volume;dx*dy];
            mass    = [mass; dx*dy*rho];
            coord   = [coord;x];
          end
          %         elseif (x(1) <= 16e3) && (x(1) >= 12e3)
          %           volume  = [volume;dx*dy];
          %           mass    = [mass; dx*dy*rho];
          %           coord   = [coord;x];
        end
      end
    end
  end
end

pCount = length(volume);

coord(:,1) = coord(:,1) - mesh.deltax*ghostCell;
coord(:,2) = coord(:,2) + mesh.deltay*ghostCell;

% stored in body1 structure

body1.volume  = volume;
body1.volume0 = volume;
body1.mass    = mass;
body1.coord   = coord;
body1.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress  = zeros(pCount,3);                % stress
body1.strain  = zeros(pCount,3);                % strain
body1.velo    = zeros(pCount,2);                % velocity
body1.deform0 = body1.deform;
C = elasticityMatrix(E,nu,'PLANE_STRAIN');
body1.C = C;
body1.gravity = g;

%% rigid body

volume=[]; mass=[]; coord=[];
for e=1:elemCount                 % start of element loop
  sctr = element(e,:);          %  element scatter vector
  pts  = node(sctr,:);
  x1 = pts(1,:); % first corner of the cell
  for i=1:ppc
    for j=1:ppc
      x(1) = x1(1) + dx*0.5 + (j-1)*dx;
      x(2) = x1(2) + dy*0.5 + (i-1)*dy;
      if (x(1) <= 16e3) && (x(1) >= 12e3) && ( x(2) >= 10e3 ) && ( x(2) <= 10e3+mesh.deltay ) 
        volume  = [volume;dx*dy];
        mass    = [mass; dx*dy*rho];
        coord   = [coord;x];
      end
    end
  end
end

pCount = length(volume);

coord(:,1) = coord(:,1) - mesh.deltax*ghostCell;
coord(:,2) = coord(:,2) + mesh.deltay*ghostCell;

% stored in body2 structure

body2.volume  = volume;
body2.volume0 = volume;
body2.mass    = mass;
body2.coord   = coord;
body2.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body2.stress  = zeros(pCount,3);                % stress
body2.strain  = zeros(pCount,3);                % strain
body2.velo    = zeros(pCount,2);                % velocity
body2.velo(:,2) = -v0;
body2.deform0 = body1.deform;
C = elasticityMatrix(E,nu,'PLANE_STRAIN');
body2.C = C;
body2.gravity = g;

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



lpx = mesh.deltax/ppc;
lpy = mesh.deltay/ppc;

%% boundary nodes

bottomNodes  = find(abs(node(:,2)-ghostCell*mesh.deltay)<1e-12);
rightNodes   = find(abs(node(:,1)-22e3-ghostCell*mesh.deltax)<1e-8);

bottomNodesGIMP = find(abs(node(:,2))<1e-12);
rightNodesGIMP  = find(abs(node(:,1)-22e3-2*ghostCell*mesh.deltax)<1e-8);

%% plot mesh, particles

figure
set(gca,'FontSize',14)
hold on
plot_mesh(node,element,'Q4','k-',1.6);
xp1 = bodies{1}.coord;
xp2 = bodies{2}.coord;
plot(xp1(:,1),xp1(:,2),'k.','markersize',15);
plot(xp2(:,1),xp2(:,2),'r.','markersize',15);
plot(node(bottomNodes,1),node(bottomNodes,2),'r.','markersize',15);
plot(node(rightNodes,1),node(rightNodes,2),'r.','markersize',15);
plot(node(bottomNodesGIMP,1),node(bottomNodesGIMP,2),'r.','markersize',15);
plot(node(rightNodesGIMP,1),node(rightNodesGIMP,2),'r.','markersize',15);

%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,2);  % nodal momentum vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector (no need for this exam)
nvelo     = zeros(nodeCount,2);  % nodal velocity vector (time t+dt)
nvelo0    = zeros(nodeCount,2);  % nodal velocity vector (time t)
nacce     = zeros(nodeCount,2);  % nodal acceleration vector

%% initial stress field
% dynamic relaxation to obtain equilibriuum

vs    = sqrt(E/rho);
dtc   = mesh.deltax/vs;

alpha=0.85;
tau   = 1e3; % c Delta t: relaxation time
dtime = 0.001;
time  = 1;
t     = 0;
tol=1e-8;

disp([num2str(toc),'   INITIAL STRESS BY DYNAMIC RELAXATION '])

gi = 0;
dg = g/1;

while ( gi < g )
  disp(['time step ',num2str(t)])
  bodies{1}.gravity=gi;
  %reset grid data
  nmass(:)     = 0;
  nmomentum(:) = 0;
  niforce(:)   = 0;
  neforce(:)   = 0;
  for ib=1:1                %% loop over bodies
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints;
    for ie=1:length(elems)         % loop over computational cells or elements
      e     = elems(ie);
      esctr = gimpElement{e};    % GIMP (extended) element connectivity
      %esctr = element(e,:);      % MPM
      %enode = node(esctr,:);     % element node coords
      mpts  = mpoints{e};        % particles inside element e
      for p=1:length(mpts)       % loop over particles
        pid    = mpts(p);
        xp     = body.coord(pid,:);
        stress = body.stress(pid,:);
        Mp     = body.mass(pid);
        vp     = body.velo(pid,:);
        Vp     = body.volume(pid);
        for i=1:length(esctr)  % loop over nodes of particle p
          id    = esctr(i);
          x     = xp - node(id,:);
          [N,dNdx]=getGIMP2D(x,mesh.deltax,mesh.deltay,lpx,lpy);
          %[N,dNdx]=getMPM2D(x,mesh.deltax,mesh.deltay);
          dNIdx = dNdx(1);
          dNIdy = dNdx(2);
          nmass(id)       = nmass(id)       + N*Mp;
          nmomentum(id,:) = nmomentum(id,:) + N*Mp*vp;
          niforce(id,1)   = niforce(id,1) - Vp*(stress(1)*dNIdx + stress(3)*dNIdy);
          niforce(id,2)   = niforce(id,2) - Vp*(stress(3)*dNIdx + stress(2)*dNIdy);
          neforce(id,2)   = neforce(id,2) - Mp*N*g;
        end
      end
    end
  end
  
  nforce = niforce + neforce;
  
  activeNodes=[bodies{1}.nodes];
  massInv = 1./nmass(activeNodes);
  % old velocity
  nvelo0(activeNodes,1) = nmomentum(activeNodes,1).*massInv;
  nvelo0(activeNodes,2) = nmomentum(activeNodes,2).*massInv;
  
  % update nodal momenta
  
  
  % dynamic relaxation using local non-viscous damping (Coetzee)
  %   for i=1:length(activeNodes)
  %     ii = activeNodes(i);
  %     fI = nforce(ii,:);
  %     vI = nvelo0 (ii,:);
  %     if (vI(1)*fI(1) <=0)
  %       alpha=0.;
  %     else
  %       alpha = 0.85;
  %     end;
  %     nforce(ii,1) = nforce(ii,1) - alpha*fI(1);
  %     if (vI(2)*fI(2) <=0)
  %       alpha=0.;
  %     else
  %       alpha = 0.85;
  %     end;
  %     nforce(ii,2) = nforce(ii,2) - alpha*fI(2);
  %     %nforce(ii,2) = nforce(ii,2) - alpha*abs(fI(2))*sign(vI(2));
  %   end
  
  nmomentum(activeNodes,:) = nmomentum(activeNodes,:) + nforce(activeNodes,:)*dtime;
  
  % dynamic relaxation with viscous damping (Sulsky)
  nmomentum(activeNodes,1) = nmomentum(activeNodes,1)./(1+tau*nmass(activeNodes));
  nmomentum(activeNodes,2) = nmomentum(activeNodes,2)./(1+tau*nmass(activeNodes));
  
  % compute updated velocities and accelerations
  for i=1:length(activeNodes)
    id = activeNodes(i);
    if nmass(id) > tol
      %nmass(id)
      nvelo(id,:) = nmomentum(id,:)/nmass(id);
      nacce(id,:) = nforce(id,:)/nmass(id);
    end
  end
  
  % boundary conditions here
  nvelo(bottomNodes,2) = 0;
  nacce(bottomNodes,2) = 0;
  nvelo(rightNodes,1)  = 0;
  nacce(rightNodes,1)  = 0;
  
  nvelo(bottomNodesGIMP,2) = 0;
  nacce(bottomNodesGIMP,2) = 0;
  nvelo(rightNodesGIMP,1)  = 0;
  nacce(rightNodesGIMP,1)  = 0;
  
  for ib=1:1
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints;
    % loop over computational cells or elements
    for ie=1:length(elems)
      e     = elems(ie);
      esctr = gimpElement{e};    % GIMP (extended) element connectivity
      %esctr = element(e,:);      % MPM
      mpts  = mpoints{e};        % particles inside element e
      for p=1:length(mpts)       % loop over particles
        pid  = mpts(p);
        xp   = body.coord(pid,:);
        Mp   = body.mass(pid);
        vp   = body.velo(pid,:);
        Vp   = body.volume(pid);
        Lp   = zeros(2,2);
        for i=1:length(esctr)
          id = esctr(i);
          vI = nvelo(id,:);
          aI = nacce(id,:);
          x     = xp - node(id,:);
          [N,dNdx]=getGIMP2D(x,mesh.deltax,mesh.deltay,lpx,lpy);
          %[N,dNdx]=getMPM2D(x,mesh.deltax,mesh.deltay);
          vp  = vp  + dtime * N*aI;
          xp  = xp  + dtime * N*vI;
          Lp  = Lp + vI'*dNdx;
        end
        
        bodies{ib}.velo(pid,:) = vp;
        bodies{ib}.coord(pid,:)= xp;
        
        % update stress last
        F       = ([1 0;0 1] + Lp*dtime)*reshape(bodies{ib}.deform(pid,:),2,2);
        bodies{ib}.deform(pid,:) = reshape(F,1,4);
        bodies{ib}.volume(pid  ) = det(F)*bodies{ib}.volume0(pid);
        dEps    = dtime * 0.5 * (Lp+Lp');
        epsTr2D = bodies{ib}.strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
        epsTr   = [epsTr2D(1) epsTr2D(2) 0 epsTr2D(3) 0 0];
        [sig,epsE,Dalg,yield_pos] = MCconstNAF(epsTr,E,nu,phi,psi,c);
        
        bodies{ib}.strain(pid,:)  = epsE([1 2 4]);
        bodies{ib}.stress(pid,:)  = sig([1 2 4]);
      end
    end
  end
  
  max(max(abs(bodies{1}.velo)))
  % update the element particle list
  for ib=1:length(1)
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
  
  % advance to the next time step
  gi = gi + dg;
end

figure(2)
xp1 = bodies{1}.coord;
scatter(xp1(:,1),xp1(:,2),85,bodies{1}.stress(:,2),'filled');
colorbar

disp([num2str(toc),'   DONE INITIAL DYNAMIC RELAXATION STEP '])

%% Slope stability starts here...



body1.stress    = bodies{1}.stress;
body1.nodes     = bodies{1}.nodes;
body1.elements  = bodies{1}.elements;
body1.mpoints   = bodies{1}.mpoints;

bodies{1} = body1;
bodyCount = length(bodies);



%% SOlver

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance
vs    = sqrt(E/rho);
dtc   = mesh.deltax/vs;
dtime = 0.001;
time  = dtime;
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
  disp(['time step ',num2str(t)])
  %reset grid data
  nmass(:)     = 0;
  nmomentum(:) = 0;
  niforce(:)   = 0;
  neforce(:)   = 0;
  for ib=1:1                       %% loop over bodies
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints;
    for ie=1:length(elems)         % loop over computational cells or elements
      e     = elems(ie);
      esctr = gimpElement{e};    % GIMP (extended) element connectivity
      %enode = node(esctr,:);     % element node coords
      mpts  = mpoints{e};        % particles inside element e
      for p=1:length(mpts)       % loop over particles
        pid    = mpts(p);
        xp     = body.coord(pid,:);
        stress = body.stress(pid,:);
        Mp     = body.mass(pid);
        vp     = body.velo(pid,:);
        Vp     = body.volume(pid);
        for i=1:length(esctr)  % loop over nodes of particle p
          id    = esctr(i);
          x     = xp - node(id,:);
          [N,dNdx]=getGIMP2D(x,mesh.deltax,mesh.deltay,lpx,lpy);
          dNIdx = dNdx(1);
          dNIdy = dNdx(2);
          nmass(id)       = nmass(id)       + N*Mp;
          nmomentum(id,:) = nmomentum(id,:) + N*Mp*vp;
          niforce(id,1)   = niforce(id,1) - Vp*(stress(1)*dNIdx + stress(3)*dNIdy);
          niforce(id,2)   = niforce(id,2) - Vp*(stress(3)*dNIdx + stress(2)*dNIdy);
          neforce(id,2)   = neforce(id,2) - Mp*N*g;
        end
      end
    end
  end
  
  
  activeNodes=[bodies{1}.nodes];
  massInv = 1./nmass(activeNodes);
  % old velocity
  nvelo0(activeNodes,1) = nmomentum(activeNodes,1).*massInv;
  nvelo0(activeNodes,2) = nmomentum(activeNodes,2).*massInv;
  
  % update nodal momenta
  nforce = niforce + neforce;
  nmomentum(activeNodes,:) = nmomentum(activeNodes,:) + nforce(activeNodes,:)*dtime;
  
  % compute updated velocities and accelerations   
  for i=1:length(activeNodes)
    id = activeNodes(i);
    if nmass(id) > tol
      %nmass(id)
      nvelo(id,:) = nmomentum(id,:)/nmass(id);
      nacce(id,:) = nforce(id,:)/nmass(id);
    end
  end
  
  % boundary conditions here
  nvelo(bottomNodes,2) = 0;
  nacce(bottomNodes,2) = 0;
  nvelo(rightNodes,1)  = 0;
  nacce(rightNodes,1)  = 0;
  
  nvelo(bottomNodesGIMP,2) = 0;
  nacce(bottomNodesGIMP,2) = 0;
  nvelo(rightNodesGIMP,1)  = 0;
  nacce(rightNodesGIMP,1)  = 0;
  
  nvelo(bodies{2}.nodes,2) = -v0;
  nvelo(bodies{2}.nodes,1) = 0;
  
  for ib=1:1
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints;
    % loop over computational cells or elements
    for ie=1:length(elems)
      e     = elems(ie);
      esctr = gimpElement{e};    % GIMP (extended) element connectivity
      enode = node(esctr,:);     % element node coords
      mpts  = mpoints{e};        % particles inside element e
      for p=1:length(mpts)       % loop over particles
        pid  = mpts(p);
        xp   = body.coord(pid,:);
        Mp   = body.mass(pid);
        vp   = body.velo(pid,:);
        Vp   = body.volume(pid);
        Lp   = zeros(2,2);
        for i=1:length(esctr)
          id = esctr(i);
          vI = nvelo(id,:);
          aI = nacce(id,:);
          x     = xp - node(id,:);
          [N,dNdx]=getGIMP2D(x,mesh.deltax,mesh.deltay,lpx,lpy);
          vp  = vp  + dtime * N*aI;
          xp  = xp  + dtime * N*vI;
          Lp  = Lp + vI'*dNdx;
        end
        
        bodies{ib}.velo(pid,:) = vp;
        bodies{ib}.coord(pid,:)= xp;
        
        % update stress last
        F       = ([1 0;0 1] + Lp*dtime)*reshape(bodies{ib}.deform(pid,:),2,2);
        bodies{ib}.deform(pid,:) = reshape(F,1,4);
        bodies{ib}.volume(pid  ) = det(F)*bodies{ib}.volume0(pid);
        dEps    = dtime * 0.5 * (Lp+Lp');
        epsTr2D = bodies{ib}.strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
        epsTr   = [epsTr2D(1) epsTr2D(2) 0 epsTr2D(3) 0 0];
        [sig,epsE,Dalg,yield_pos] = MCconstNAF(epsTr,E,nu,phi,psi,c);
        
        bodies{ib}.strain(pid,:)  = epsE([1 2 4]);
        bodies{ib}.stress(pid,:)  = sig([1 2 4]);
      end
    end
  end
  
  % update the rigid body
  bodies{2}.coord = bodies{2}.coord + dtime* bodies{2}.velo;
  
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
  
  %   ta = [ta;t];
  %   ka = [ka;k];
  %   sa = [sa;u];
  
  % VTK output
  
  if (  mod(istep-1,interval) == 0 )
    xp = [bodies{1}.coord];
    s  = [bodies{1}.stress];
    stress = [s sum(s,2)/3];
    data.stress = stress;
    data.velo=[bodies{1}.velo];
    vtkFile = sprintf('%s%d',vtkFileName,istep-1);
    VTKParticles(xp,vtkFile,data);
  end
  
  % advance to the next time step
  
  t     = t + dtime;
  istep = istep + 1;
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

disp([num2str(toc),'   DONE '])
