% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Particles are generated as either Gauss points of the background mesh or
% as points that regularly divide the grid cells.
%
% New format: particles are grouped in bodies to facilitate multibody
% simulations. Do not loop over grid cells, instead, loop directly over the
% particles.
%
% Some parts are done using MEX files in an attempt to speed-up the code.
% This serves a good example to learn programming MEX files.
%
%
% Two elastic disks come into contact.
% Example taken from Sulsky et al. 1994 paper.
%
% Vinh Phu Nguyen
% 25 June 2014, Saigon, Vietnam.

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
addpath ../../mex/

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

fac   = 1; % 1: explicit time
% 100: implicit is better
v     = 0.1/fac;   % initial particle velocity

%useMex = 1; % use MEX functions (done in 21seconds)
useMex = 0; % do not use MEX functions (done in 314s)

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell= 0;
lx       = 1;
ly       = 1;
numx2    = 30;       % number of elements along X direction
numy2    = 30;      % number of elements along Y direction
[mesh]   = buildGrid2D(lx,ly,numx2,numy2, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;

%% particle generation
% body1

ppc    = 3; % # of particle per cell is ppc x ppc
useGPs = 1; % particles at Gauss points
useGPs = 1; % particles at geometric points regularly divide the elements

if (useGPs==0)
  dx = mesh.deltax/(ppc);
  dy = mesh.deltay/(ppc);
end

[W,Q]=quadrature(  3, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

center = [0.2 0.2];
radius = 0.2;

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
      r  = norm(x-center);
      if ( r-radius < 0 )
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
        r  = norm(x-center);
        if ( r-radius < 0 )
          volume  = [volume;dx*dy];
          mass    = [mass; dx*dy*rho];
          coord   = [coord;x];
        end
      end
    end
  end
end

pCount = length(volume);

% stored in body1 structure

body1.volume  = volume;
body1.volume0 = volume;
body1.mass    = mass;
body1.coord   = coord;
body1.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress  = zeros(pCount,3);                % stress
body1.strain  = zeros(pCount,3);                % strain
body1.velo    = ones(pCount,2)*v;               % velocity
body1.C       = C;
body1.deform0 = body1.deform;

% body 2
center = [0.8 0.8];
radius = 0.2;

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
      r  = norm(x-center);
      if ( r-radius < 0 )
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
        r  = norm(x-center);
        if ( r-radius < 0 )
          volume  = [volume;dx*dy];
          mass    = [mass; dx*dy*rho];
          coord   = [coord;x];
        end
      end
    end
  end
end

body2.volume  = volume;
body2.volume0 = volume;
body2.mass    = mass;
body2.coord   = coord;
body2.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body2.stress  = zeros(pCount,3);                % stress
body2.strain  = zeros(pCount,3);                % strain
body2.velo    = -ones(pCount,2)*v;              % velocity
body2.C       = C;
body2.deform0 = body2.deform;

bodies    = cell(2,1);
bodies{1} = body1;
bodies{2} = body2;
bodyCount = length(bodies);

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

bodies = findActiveElemsAndNodes(bodies,mesh);

%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,2);  % nodal momentum vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector (no need for this exam)
nvelo = zeros(nodeCount,2);      % nodal velocity vector
nacce = zeros(nodeCount,2);      % nodal acceleration vector

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

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.3*mesh.deltax/c;
%time  = 60; % low velocity, used to test implicit
time  = 2;% explicit, high velocity
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
  disp(['time step ',num2str(t)])
  if (useMex==0)
    %reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    for ib=1:bodyCount                         %% loop over bodies
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
        end
      end
    end
  else
    % MEX function
    [nmass,nmomentum,niforce] = ParticlesToNodes(bodies,mesh);
  end
  
  % update nodal momenta
  
  activeNodes=[bodies{1}.nodes; bodies{2}.nodes];
  nmomentum(activeNodes,:) = nmomentum(activeNodes,:) + niforce(activeNodes,:)*dtime;
  
  % compute updated velocities and accelerations
  nvelo(activeNodes,1) = nmomentum(activeNodes,1)./nmass(activeNodes);
  nvelo(activeNodes,2) = nmomentum(activeNodes,2)./nmass(activeNodes);
  
  nacce(activeNodes,1) = niforce(activeNodes,1)./nmass(activeNodes);
  nacce(activeNodes,2) = niforce(activeNodes,2)./nmass(activeNodes);
  % note that there is no boundary conditions for this example
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
        [stress,strain] = updateStressHooke(dStrain,C,bodies{ib}.stress(p,:),bodies{ib}.strain(p,:));
        
        bodies{ib}.strain(p,:)  = strain;
        bodies{ib}.stress(p,:)  = stress;
      end
    end
  else
    UpdateParticles(bodies,mesh,nvelo,nacce,dtime); % MEX function
  end
  
  bodies{1}.coord0 = bodies{1}.coord;
  bodies{2}.coord0 = bodies{2}.coord;
  
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
  
  bodies = findActiveElemsAndNodes(bodies,mesh);
  
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

%ss=load('mpmTwoDisksMex.mat');

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(ta(1:end),sa(1:end),'r-','LineWidth',2);
plot(ss.ta(1:end),ss.ka(1:end),'b.','LineWidth',1.6);
plot(ss.ta(1:end),ss.sa(1:end),'r.','LineWidth',2);
%plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy (x1E-3)')
legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%axis([0 3 0 3])

disp([num2str(toc),'   DONE '])
