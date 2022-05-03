% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Particles are generated using an unstructured FE mesh (gmsh).
%
% Shape functions: using standard FE shape function routine by converting
% particle position to natural coordinates.
%
% Implicit time integration: Newmark scheme.
%
% Two elastic disks come into contact.
% Example taken from Sulsky et al. 1994 paper.
%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia.
% 10 August 2015.

%%

addpath ../../grid/
addpath ../../basis/
addpath ../../particleGen/
addpath ../../constitutiveModels/
addpath ../../util/
addpath ../../geoMesh/
addpath ../../externals/
addpath ../../postProcessing/

%%
clc
clear
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Define some constants: material properties, time integration constants...
%
E   = 1000;        % Young's modulus
nu  = 0.3;         % Poisson ratio
rho = 1000;        % density
kappa = 3-4*nu;    % Kolosov constant
mu    = E/2/(1+nu);% shear modulus
v0    = 0.1;       % initial particle velocity
bulk  = E/3/(1-2*nu);

% Newmark parameters
beta  = 0.25;
gamma = 0.5;

interval     = 10;% time interval for saving vtp files.
vtkFileName  = 'mpm2DTwoDisks';

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell= 0;
lx       = 1;
ly       = 1;
numx2    = 20;       % number of elements along X direction
numy2    = 20;      % number of elements along Y direction
[mesh]   = buildGrid2D(lx,ly,numx2,numy2, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;

%%
ppc           = [2 2];
cir1.center   = [0.2 0.2];
cir1.radius   = 0.2;
cir2.center   = [0.8 0.8];
cir2.radius   = 0.2;
[res1]        = generateMPForCircle(cir1,ppc,mesh);
[res2]        = generateMPForCircle(cir2,ppc,mesh);

pCount  = 2*size(res1.position,1);
volume  = [res1.volume;res2.volume];
volume0 = volume;
mass    = volume*rho;
coords  = [res1.position;res2.position];
coords0 = coords;
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stress
strain  = zeros(pCount,3);                % strain
velo    = zeros(pCount,2);                % velocity
density = rho*ones(pCount,1);             % density

% initial velocities, initial stress=0
for p=1:pCount
  if coords(p,1) < 0.5
    velo(p,:) = [v0 v0];
  else
    velo(p,:) = [-v0 -v0];        
  end
end

dis     = zeros(pCount,2);                    % displacement
dis0    = zeros(pCount,2);                    % old displacement
velo0   = velo;                               % old velocity
acce    = zeros(pCount,2);                    % acceleration
acce0   = zeros(pCount,2);                    % old acceleration
%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

pElems  = ones(pCount,1);
mpoints = cell(elemCount,1);

for p=1:pCount
  x = coords(p,1);
  y = coords(p,2);
  e = floor(x/mesh.deltax) + 1 + numx2*floor(y/mesh.deltay);
  pElems(p) = e;
end

for e=1:elemCount
  id  = find(pElems==e);
  mpoints{e}=id;
end

%% node quantities

dofCount  = 2*nodeCount;

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,2);  % nodal momentum vector
nmomentum1= zeros(nodeCount,2);  % nodal momentum vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector (no need for this exam)
% note that nodal velo/disp/acce stored as vectors not matrices
nvelo     = zeros(dofCount,1);   % nodal velocity
ndis      = zeros(dofCount,1);   % nodal displacement
nacce     = zeros(dofCount,1);   % nodal acceleration

ndis0     = zeros(dofCount,1);   % old nodal displacement
nvelo0    = zeros(dofCount,1);   % old nodal velocity
nacce0    = zeros(dofCount,1);   % old nodal acceleration

activeNodes = unique(element(unique(pElems),:));
activeDofs(1:2:2*length(activeNodes)) = activeNodes*2-1;
activeDofs(2:2:2*length(activeNodes)) = activeNodes*2;

K         = sparse(dofCount,dofCount); % stiffness matrix
M         = sparse(dofCount,dofCount); % mass matrix

%% plot mesh, particles

figure(1)
hold on
%plot_mesh(node1,element1,elemType,'r-',0.2);
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
axis off

% save some quantities for post-processing

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.02;
time  = 2.5;%10*dtime;
t     = 0;
dtime2    = dtime*dtime;
fac       = 1/(beta*dtime2);

nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
  disp(['time step ',num2str(t)])
  % reset grid data
  nmass(:)     = 0;
  ndis0(:)     = 0;
  nacce0(:)    = 0;
  nvelo0(:)    = 0;
  nmomentum(:) = 0; nmomentum1(:) = 0;
  niforce(:)   = 0;
  K(:,:)       = 0; M(:,:)       = 0;
  %     dis(:)       = 0;
  %     acce(:)      = 0;
  % loop over computational cells or elements
  for e=1:elemCount
    esctr = element(e,:);      % element connectivity
    enode = node(esctr,:);     % element node coords
    mpts  = mpoints{e};        % particles inside element 'e'
    ncount= length(esctr);     % number of nodes of element 'e'
    for p=1:length(mpts)       % loop over particles
      pid    = mpts(p);        % index of particle 'p'
      sigma  = stress(pid,:);
      vol    = volume(pid  );
      m      = mass  (pid  );
      vel    = velo  (pid,:);
      d      = dis0  (pid,:);
      acc    = acce0 (pid,:);
      xpa    = coords(pid,:);
      % evaluate basis functions and derivatives of 4 nodes at
      % particle 'p'
      pt(1)= (2*xpa(1)-(enode(1,1)+enode(2,1)))/mesh.deltax;
      pt(2)= (2*xpa(2)-(enode(2,2)+enode(3,2)))/mesh.deltay;
      
      [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
      J0       = enode'*dNdxi;             % element Jacobian matrix
      invJ0    = inv(J0);
      dNdx     = dNdxi*invJ0;
      % particle to node
      for i=1:ncount
        id    = esctr(i);                 % index of node 'i'
        BI    = [dNdx(i,1) 0;0 dNdx(i,2);dNdx(i,2) dNdx(i,1)];
        nmass(id)       = nmass(id)        + N(i)*m;
        nmomentum(id,:) = nmomentum(id,:)  + N(i)*m*vel;
        nmomentum1(id,:)= nmomentum1(id,:) + N(i)*m*acc;
        %niforce(id,:)   = niforce(id,:)    - vol*sigma*BI;
        %ndis0([2*id-1;2*id])  = ndis0([2*id-1;2*id])   + N(i)*d';
        %nacce0([2*id-1;2*id]) = nacce0([2*id-1;2*id])  + N(i)*acc';
        for j=1:ncount
          jd    = esctr(j);
          BJ    = [dNdx(j,1) 0;0 dNdx(j,2);dNdx(j,2) dNdx(j,1)];
          K([2*id-1 2*id],[2*jd-1 2*jd]) = K([2*id-1 2*id],[2*jd-1 2*jd]) + ...
            vol*BI'*C*BJ;
          
          M([2*id-1 2*id],[2*jd-1 2*jd]) = M([2*id-1 2*id],[2*jd-1 2*jd]) + ...
            [m*N(i)*N(j) 0;0 m*N(i)*N(j) ];
        end
      end
    end
  end
  
  nvelo0(2*activeNodes-1) = nmomentum(activeNodes,1)./nmass(activeNodes);
  nvelo0(2*activeNodes  ) = nmomentum(activeNodes,2)./nmass(activeNodes);
  % aI = sum_p m_p a_p /m_I
  nacce0(2*activeNodes-1) = nmomentum1(activeNodes,1)./nmass(activeNodes);
  nacce0(2*activeNodes  ) = nmomentum1(activeNodes,2)./nmass(activeNodes);
  
  %% update nodal displacement by solving Ax=b
  % ndisp0 = 0: as grid never moves.  
  dtil      = dtime*nvelo0 + 0.5*dtime2*(1-2*beta)*nacce0;
  A         = fac*M + K;
  f         = fac*M*dtil;
  
  % apply boundary conditions
  freeNodes= setdiff(1:nodeCount,activeNodes)';
  udofs    = 2*freeNodes-1;
  vdofs    = 2*freeNodes;
  uFixed   = zeros(1,length(freeNodes));
  vFixed   = zeros(1,length(vdofs));
  [A,f]    = applyDirichletBCs(A,f,udofs,vdofs,uFixed,vFixed);
  % solve the system to get updated nodal displacement
  ndis     = A\f;
  
  % update nodal acceleration and velocity
  nacce = (4/dtime2)*ndis - (4/dtime)*nvelo0 - nacce0;
  nvelo = nvelo0 + (1-gamma)*dtime*nacce0 + gamma*dtime*nacce;
  
  %% update particle displacement and acceleration and stresses
  k = 0; u = 0;
  for e=1:elemCount
    esctr = element(e,:);
    enode = node(esctr,:);
    mpts  = mpoints{e};
    ncount= length(esctr);
    for p=1:length(mpts) % loop over particles
      pid  = mpts(p);
      xpa  = coords(pid,:);
      Lp   = zeros(2,2);
      d    = [0 0];
      a    = [0 0];
      for i=1:ncount
        id           = esctr(i);
        idx          = [2*id-1;2*id];
        vI           = nvelo(idx);
        aI           = nacce(idx);
        dI           = ndis (idx);
        x            = xpa - node(id,:);
        [N,dNdx]     = getMPM2D(x,mesh.deltax,mesh.deltay);
        d            = d   +  N*dI';
        a            = a   +  N*aI';
        Lp           = Lp  +  vI*dNdx;
      end
      
      dis (pid,:)  = d;
      acce(pid,:)  = a;
      
      % update particle position and velocity
      velo(pid,:)   = velo0(pid,:)   +  gamma*dtime*a + (1-gamma)*dtime*acce0(pid,:);
      coords(pid,:) = coords(pid,:)  +  d;
      
      k    = k + 0.5*(velo(pid,1)^2+velo(pid,2)^2)*mass(pid);
      u    = u + 0.5*volume(pid)*stress(pid,:)*strain(pid,:)';
      
      F              = ([1 0;0 1] + Lp*dtime)*reshape(deform(pid,:),2,2);
      deform(pid,:)  = reshape(F,1,4);
      volume(pid)    = det(F)*volume0(pid);
      dEps           = dtime * 0.5 * (Lp+Lp'); % strain matrix
      dsigma         = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
      stress(pid,:)  = stress(pid,:) + dsigma';
      strain(pid,:)  = strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
    end
  end
  
  dis0  = dis;
  acce0 = acce;
  velo0 = velo;
  
  % store time,velocty for plotting
  %
  %     pos{istep} = coords;
  %     vel{istep} = velo;
  
  % update the element particle list
  
  for p=1:pCount
    x = coords(p,1);
    y = coords(p,2);
    e = floor(x/mesh.deltax) + 1 + numx2*floor(y/mesh.deltay);
    pElems(p) = e;
  end
  
  for e=1:elemCount
    id  = find(pElems==e);
    mpoints{e}=id;
  end
  
  activeNodes         = unique(element(unique(pElems),:));
  activeDofs=[];
  activeDofs(1:2:2*length(activeNodes)) = activeNodes*2-1;
  activeDofs(2:2:2*length(activeNodes)) = activeNodes*2;
  
  % store time,velocty for plotting
  
  ta = [ta;t];
  ka = [ka;k];
  sa = [sa;u];
  
  % VTK output
  
%   if (  mod(istep,interval) == 0 )
%     xp = coords;
%     vtkFile = sprintf('../results/%s%d',vtkFileName,istep);
%     data.stress  = [stress zeros(size(stress,1),1)];
%     data.velo    = velo;
%     VTKParticles(xp,vtkFile,data);
%   end
  
  
  % advance to the next time step
  
  t     = t + dtime;
  istep = istep + 1;
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

ss=load('mpm-explicit-2disks-dt001.mat');
%

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'red-','LineWidth',1.6);
%plot(ta(1:end),sa(1:end),'r-','LineWidth',2);
%plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
plot(ss.ta(1:end),ss.ka(1:end),'black-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy')
legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
axis([0 3.5 0 3])

%savefile = 'mpm-implicit-2disks-dt002.mat';
%save(savefile,'ta','ka','sa');

disp([num2str(toc),'   DONE '])

%% plot curves for different time steps

% ss1 = load('mpm-implicit-2disks-dt001.mat');
% ss2 = load('mpm-implicit-2disks-dt002.mat');
%
% figure
% set(gca,'FontSize',14)
% hold on
% plot(ss1.ta(1:end),ss1.ka(1:end),'black-','LineWidth',2.1);
% plot(ss1.ta(1:end),ss1.sa(1:end),'red-','LineWidth',2.1);
% plot(ss1.ta(1:end),ss1.ka(1:end)+ss1.sa(1:end),'g-','LineWidth',2.1);
% xlabel('Time')
% ylabel('Energy')
% legend('kinetic','strain','total')

%exportfig(gcf,'implicit-sulsky-2disks-dt01.eps',opts)

% write the background gid

% Ux= zeros(size(node,1),1);
% Uy= zeros(size(node,1),1);
% sigmaXX = zeros(size(node,1),1);
% sigmaYY = zeros(size(node,1),1);
% sigmaXY = zeros(size(node,1),1);
%
% VTKPostProcess(node,element,2,'Quad4','grid',...
%     [sigmaXX sigmaYY sigmaXY],[Ux Uy]);


