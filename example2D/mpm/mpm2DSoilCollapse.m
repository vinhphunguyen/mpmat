% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Particles are generated as either Gauss points of the background mesh or
% as points that regularly divide the grid cells.
%
% Collapse of a soild column with frictional boundary conditions.
% Soil: Morh-Coulomb model.
%
% 10 August: node explosion without double mapping.
%
% Vinh Phu Nguyen
% 8 August 2014, Adelaide, Australia.

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

vtkFileName  = 'mpm2DSoilCollapse';
interval     = 50;

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell= 0;
lx       = 24e3;
ly       = 8e3;
numx2    = 60;       % number of elements along X direction
numy2    = 20;      % number of elements along Y direction
[mesh]   = buildGrid2D(lx,ly,numx2,numy2, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;

%% particle generation
% body1

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
        if ( x(1) <= 4e3 )
          volume  = [volume;dx*dy];
          mass    = [mass; dx*dy*rho];
          coord   = [coord;x];
        end
      end
    end
  end
end

pCount = length(volume);

coord(:,1) = coord(:,1) + 10e3;

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

% initial stress field

for p=1:length(body1.mass)
  y = body1.coord(p,2);
  h = 8e3 - y;
  sigmayy = -g*rho*h;
  sigmaxx = K0*sigmayy;
  body1.stress(p,[1 2])=[sigmaxx sigmayy];
end

bodies    = cell(1,1);
bodies{1} = body1;
bodyCount = length(bodies);

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

bodies = findActiveElemsAndNodes(bodies,mesh);

%% boundary nodes

bottomNodes = find(node(:,2)<1e-12);

%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,2);  % nodal momentum vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector (no need for this exam)
nvelo     = zeros(nodeCount,2);  % nodal velocity vector (time t+dt)
nvelo0    = zeros(nodeCount,2);  % nodal velocity vector (time t)
nacce     = zeros(nodeCount,2);  % nodal acceleration vector

%% plot mesh, particles

figure
set(gca,'FontSize',14)
hold on
plot_mesh(node,element,'Q4','k-',1.6);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','cy-',1.6);
xp1 = bodies{1}.coord;
plot(xp1(:,1),xp1(:,2),'k.','markersize',15);

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

dtime = 0.001;
time  = 0.98;
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
  disp(['time step ',num2str(t)])
  if (1)
    %reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    neforce(:)   = 0;
    for ib=1:bodyCount                %% loop over bodies
      body      = bodies{ib};
      for p=1:length(body.mass)       % loop over particles
        xp     = body.coord(p,:);
        stress = body.stress(p,:);
        Mp     = body.mass(p);
        vp     = body.velo(p,:);
        Vp     = body.volume(p);
        nodes  = getNodesForParticle2D(xp(1), xp(2), mesh.deltax, mesh.deltay, mesh.numx);
        for i=1:length(nodes)  % loop over nodes of particle p
          id    = nodes(i);
          x     = xp - node(id,:);
          [N,dNdx]=getMPM2D(x,mesh.deltax,mesh.deltay);
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
  else
    % MEX function
    [nmass,nmomentum,niforce] = ParticlesToNodes(bodies,mesh);
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
  nvelo(activeNodes,1) = nmomentum(activeNodes,1).*massInv;
  nvelo(activeNodes,2) = nmomentum(activeNodes,2).*massInv;
  
  nacce(activeNodes,1) = nforce(activeNodes,1).*massInv;
  nacce(activeNodes,2) = nforce(activeNodes,2).*massInv;
  
  % frictional boundary conditions here
  
  nvelo(bottomNodes,2) = 0;
  nacce(bottomNodes,2) = 0;
  
  contactNodes = intersect(bottomNodes,activeNodes);
  
  for in=1:length(contactNodes)
    id       =  contactNodes(in);
    velo1    = nvelo(id,:);
    
    velocm   = [0 0];
    
    nI       = [0 -1];    
    deltaVe  = velo1 - velocm;
    D        = dot(deltaVe, nI);
    C        = abs(deltaVe(1)*nI(2) - deltaVe(2)*nI(1));
    muPrime  = min(mu,C/D);
    if ( D >= 0 )
      nvelo(id,:) = velo1 - D*(nI+(1/C)*muPrime*[nI(2)*(deltaVe(1)*nI(2) - deltaVe(2)*nI(1)) -nI(1)*(deltaVe(1)*nI(2) - deltaVe(2)*nI(1))]);
      nacce(id,:) = (1/dtime)*( nvelo(id,:) - nvelo0(id,:) );
      % disp('approaching')
    else
      % disp('separating')
    end
  end
  
  if (useMex==0)
    for ib=1:bodyCount
      body      = bodies{ib};
      for p=1:length(body.mass)       % loop over particles
        xp   = body.coord(p,:);
        xp0  = body.coord(p,:);
        vp   = body.velo(p,:);
        nodes  = getNodesForParticle2D(xp(1), xp(2), mesh.deltax, mesh.deltay, mesh.numx);
        Lp   = zeros(4,1);
        for i=1:length(nodes)
          id  = nodes(i);
          x   = xp0 - node(id,:);
          [N,dNdx]=getMPM2D(x,mesh.deltax,mesh.deltay);
          
          xp  = xp + dtime * N*nvelo(id,:);%nmomentum(id,:)/nmass(id);
          vp  = vp + dtime * N*nacce(id,:);%niforce(id,:)/nmass(id);
          Lp(1) = Lp(1) + nvelo(id,1)*dNdx(1);
          Lp(2) = Lp(2) + nvelo(id,1)*dNdx(2);
          Lp(3) = Lp(3) + nvelo(id,2)*dNdx(1);
          Lp(4) = Lp(4) + nvelo(id,2)*dNdx(2);
        end
        bodies{ib}.coord(p,:)= xp;
        bodies{ib}.velo(p,:) = vp;
        
        a11 = 1.+dtime*Lp(1); a12 = dtime*Lp(2); a21 = dtime*Lp(3); a22 = 1.+dtime*Lp(4);
        fxx = body.deform(p,1); fyx = body.deform(p,2);
        fxy = body.deform(p,3); fyy = body.deform(p,4);
        Fxx = a11*fxx + a12*fyx;
        Fxy = a11*fxy + a12*fyy;
        Fyx = a21*fxx + a22*fyx;
        Fyy = a21*fxy + a22*fyy;
        bodies{ib}.deform(p,1) = Fxx; bodies{ib}.deform(p,2) = Fyx;
        bodies{ib}.deform(p,3) = Fxy; bodies{ib}.deform(p,4) = Fyy;
        detF     = Fxx*Fyy - Fxy*Fyx;
        bodies{ib}.volume(p) = bodies{ib}.volume0(p) * detF;
        dStrain = dtime*[Lp(1) Lp(4) Lp(2)+Lp(3)]; % Voigt notation
        epsTr2D = [bodies{ib}.strain(p,:)+dStrain];
        epsTr   = [epsTr2D(1) epsTr2D(2) 0 epsTr2D(3) 0 0];
        [sig,epsE,Dalg,yield_pos] = MCconstNAF(epsTr,E,nu,phi,psi,c);
        
        bodies{ib}.strain(p,:)  = epsE([1 2 4]);
        bodies{ib}.stress(p,:)  = sig([1 2 4]);
      end
    end
  else
    UpdateParticles(bodies,mesh,nvelo,nacce,dtime); % MEX function
  end
  
  % update the element particle list
  
  bodies = findActiveElemsAndNodes(bodies,mesh);
  
  % store time,velocty for plotting
  
  %   ta = [ta;t];
  %   ka = [ka;k];
  %   sa = [sa;u];
  
  % advance to the next time step
  
  t     = t + dtime;
  istep = istep + 1;
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

disp([num2str(toc),'   DONE '])
