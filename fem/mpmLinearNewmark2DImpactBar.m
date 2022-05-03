% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Particles are generated using an unstructured FE mesh (gmsh).
%
% Shape functions: using standard FE shape function routine by converting
% particle position to natural coordinates.
%
% Implicit time integration: Newmark scheme.
% Wave propagation in 1D bar.
%
% Vinh Phu Nguyen
% South Australia.
% 2 January 2016.

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
E   = 1;        % Young's modulus
nu  = 0.;       % Poisson ratio
rho = 1;        % density
sigma0 = 1.;
c      = sqrt(E/rho);

% Newmark parameters
beta  = 0.25;
gamma = 0.5;

CFL   = 0.5;

stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell= 0;
lx       = 10;
ly       = 1;

numx2    = 100;       % number of elements along X direction
numy2    = 1;         % number of elements along Y direction
[mesh]   = buildGrid2D(lx,ly,numx2,numy2, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;

%% Boundary nodes (for essential BCs)

fixedNodes  = mesh.rNodes; % nodes
forcedNodes = mesh.lNodes; % nodes
uFixed     = zeros(1,length(fixedNodes));
vFixed     = zeros(1,length(fixedNodes));
udofs      = 2*fixedNodes-1;
vdofs      = 2*fixedNodes;

%%
ppc           = [2 2];
rec.x         = [0 lx];
rec.y         = [0 ly];
[res1]        = generateMPForRectangle(rec,ppc,mesh);

[W,Q]   = quadrature( 2, 'GAUSS', 2 ); % 2x2 Gauss point rule

volume = [];
pos    = [];

for e=1:elemCount
  esctr = element(e,:);
  enode = node(esctr,:);
  % loop over Gauss point
  for p=1:length(W)
    pt       = Q(p,:);
    [N,dNdxi]= lagrange_basis('Q4',pt);   % element shape functions
    J0       = enode'*dNdxi;                  % element Jacobian matrix
    detJ     = det(J0);
    vol      = detJ * W(p);
    volume   = [volume;vol];
    pos      = [pos;N'*enode];
  end
end

res1.volume   = volume;
res1.position = pos;

pCount  = size(res1.position,1);
volume  = res1.volume;
volume0 = volume;
mass    = volume*rho;
coords  = res1.position;
coords0 = coords;
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stress
strain  = zeros(pCount,3);                % strain
velo    = zeros(pCount,2);                % velocity
density = rho*ones(pCount,1);             % density

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
nmomentum2= zeros(nodeCount,2);  % nodal momentum vector

niforce   = zeros(dofCount,1);   % nodal internal force vector
neforce   = zeros(dofCount,1);   % nodal external force vector 
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

K         = zeros(dofCount,dofCount); % stiffness matrix
M         = zeros(dofCount,dofCount); % mass matrix

%% plot mesh, particles

figure(1)
hold on
%plot_mesh(node1,element1,elemType,'r-',0.2);
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
plot(node(fixedNodes,1),node(fixedNodes,2),'bs','MarkerSize',9,'MarkerFaceColor','blue')
plot(node(forcedNodes,1),node(forcedNodes,2),'rs','MarkerSize',9,'MarkerFaceColor','red')
axis off

% save some quantities for post-processing

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

dtime  = CFL*mesh.deltax/c;
dtime2 = dtime*dtime;
time   = 0.5*lx/c;20*dtime;
t      = 0;
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
  nmomentum(:) = 0; nmomentum1(:) = 0; nmomentum2(:) = 0;
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
        nmomentum2(id,:)= nmomentum2(id,:) + N(i)*m*d;
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

  ndis0(2*activeNodes-1) = nmomentum2(activeNodes,1)./nmass(activeNodes);
  ndis0(2*activeNodes  ) = nmomentum2(activeNodes,2)./nmass(activeNodes);
  
  %% update nodal displacement by solving Ax=b
  % ndisp0 = 0: as grid never moves.
  neforce(2*forcedNodes-1) = 0.5*sigma0;
  fac       = 1/(beta*dtime2);
  dtil      = ndis0 + dtime*nvelo0 + 0.5*dtime2*(1-2*beta)*nacce0;
  A         = fac*M + K;
  f         = fac*M*dtil + neforce;
  
  % apply boundary conditions
  [A,f]    = applyDirichletBCs(A,f,udofs,vdofs,uFixed,vFixed);
  % solve the system to get updated nodal displacement
  ndis     = A\f;
  
  % update nodal acceleration and velocity
  nacce = fac*(ndis-dtil);
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
      %coords(pid,:) = coords(pid,:)  +  d;
      
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

numSig = zeros(2*mesh.elemCount,1);
numLoc = zeros(2*mesh.elemCount,1);
ii = 1;
for i=1:elemCount
% 1 => 1 2
% 2 => 5 6
% 4*(i-1)+1
numLoc(ii)   = coords(4*(i-1)+1,1);
numLoc(ii+1) = coords(4*(i-1)+2,1);

numSig(ii)   = stress(4*(i-1)+1,1);
numSig(ii+1) = stress(4*(i-1)+2,1);
ii = ii + 2;
end

% exact solution
x        = linspace(0,lx,100);
exactSig = zeros(length(x),1);
for i=1:length(x)
  exactSig(i) = -sigma0*heaviside(c*time-x(i));
end

figure
hold on
plot(numLoc/lx,numSig,'black-','LineWidth',1.8)
plot(x/lx,exactSig,'blue-','LineWidth',1.8)
xlabel('x/L')
ylabel('stress')
set(gca,'FontSize',16)
legend('FEM','Exact')
grid on


disp([num2str(toc),'   DONE '])

%% plot curves for different time steps




