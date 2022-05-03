% This file implements the  coupled IGA-Material Point Method
%
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
%
% Collision of two elastic disks.
%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 14 September 2015.

addpath ../../grid/
addpath ../../basis/
addpath ../../particleGen/
addpath ../../constitutiveModels/
addpath ../../util/
addpath ../../geoMesh/
addpath ../../externals/
addpath ../../postProcessing/
addpath ../../nurbs/nurbs-util/
addpath ../../nurbs/nurbs-geopdes/inst/

%%
clc
clear 
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

tic;
%%
E   = 1000;        % Young's modulus
nu  = 0.3;         % Poisson ratio
rho = 1000;        % density
stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C  = elasticityMatrix(E,nu,stressState);
v0 = 0.1;

%% Computational grid (prefer square cells)
l = 1.5;
w = 1.5;

noY0      = 30;         % number of elements along Y direction
noX0      = 30;         % number of elements along X direction
ghostCell = 0;

[bGrid]    = buildGrid2D(l,w,noX0,noY0, ghostCell);

bGrid.node(:,1) = bGrid.node(:,1) - 0.3;
bGrid.node(:,2) = bGrid.node(:,2) - 0.3;

node      = bGrid.node;
element   = bGrid.element;
deltax    = bGrid.deltax;
deltay    = bGrid.deltay;
elemCount = bGrid.elemCount;
nodeCount = bGrid.nodeCount;
numx2     = bGrid.numx;
numy2     = bGrid.numy;
Vc        = deltax * deltay;
xmin      = min(node(:,1));
ymin      = min(node(:,2));
%% particles as a NURBS surface


%
radius   = 0.2;
center1  = [0.2 0.2];
center2  = [0.8 0.2];
circle1  = buildNURBSCircle(radius,center1);
circle2  = buildNURBSCircle(radius,center2);
% perform refinement on NURBS surfaces
refineCount = 2;
circle1 = hRefineNURBS(circle1,refineCount);
circle2 = hRefineNURBS(circle2,refineCount);
% convert NURBS to a FEM-like format
mesh1    = generateIGA2DMesh(circle1);
mesh2    = generateIGA2DMesh(circle2);

% Gauss quadrature rule
noGPs = 2;
[W,Q] = quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

% Pre-compute Bernstein basis and derivatives for ONE Bezier element
noBasis = (mesh1.p+1)*(mesh1.q+1);
noGpEle = noGPs*noGPs;

shapes  = zeros(noGpEle,noBasis);
derivs  = zeros(noGpEle,noBasis,2);

for gp=1:size(W,1)
  [shapes(gp,:), derivs(gp,:,:)] = getShapeGradBernstein2D(mesh1.p,mesh1.q,Q(gp,1),Q(gp,2));
end


pCount         = mesh1.nodeCount;
body1.coord    = mesh1.node;
body1.coord0   = mesh1.node;
body1.deform   = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress   = zeros(pCount,3);                % stress
body1.strain   = zeros(pCount,3);                % strain
body1.velo     = zeros(pCount,2);                % velocity
body1.up       = zeros(pCount,2);                % velocity
body1.velo(:,1) = v0;
%body1.velo(:,2) = v0;
body1.fp       = zeros(pCount,2);                
body1.mass     = zeros(pCount,1);    
body1.volume   = zeros(pCount,1); 
body1.mesh     = mesh1;
body1.nurbs    = circle1;
body1.stress   = zeros(body1.mesh.elemCount,noGpEle,3);
body1.strain   = zeros(body1.mesh.elemCount,noGpEle,3);

pCount         = mesh2.nodeCount;
body2.coord    = mesh2.node;
body2.coord0   = mesh2.node;
body2.deform   = repmat([1 0 0 1],pCount,1);     % gradient deformation
body2.stress   = zeros(pCount,3);                % stress
body2.strain   = zeros(pCount,3);                % strain
body2.velo     = zeros(pCount,2);                % velocity
body2.up       = zeros(pCount,2);                % velocity
body2.velo(:,1) = -v0;
%body2.velo(:,2) = -v0;
body2.fp       = zeros(pCount,2);                
body2.mass     = zeros(pCount,1);    
body2.volume   = zeros(pCount,1);
body2.mesh     = mesh2;
body2.nurbs    = circle2;
body2.stress   = zeros(body2.mesh.elemCount,noGpEle,3);
body2.strain   = zeros(body2.mesh.elemCount,noGpEle,3);

% determine particle mass
M = zeros(mesh1.nodeCount,mesh1.nodeCount);
for e=1:mesh1.elemCount
  sctr   = mesh1.element(e,:);         %  element scatter vector
  nn     = length(sctr);
  Ce     = mesh1.bezierExtractor(:,:,e); % element Bezier extraction operator
  we     = diag(mesh1.weights(sctr));    % element weights
  pts    = mesh1.node(sctr,:);           % element nodes
  Wb     = Ce'*mesh1.weights(sctr);      % element Bezier weights    
  for gp=1:size(W,1)                    % loop over Gauss points
    pt      = Q(gp,:);
    wt      = W(gp);
    % Bernstein basis and derivatives at GP gp
    Be      = shapes(gp,:)';
    dBedxi  = reshape(derivs(gp,:,:),noBasis,2);
    % Bezier weight functions (denomenator of NURBS)
    wb        = dot(Be,Wb);
    dwbdxi(1) = dot(dBedxi(:,1),Wb);
    dwbdxi(2) = dot(dBedxi(:,2),Wb);
    % Shape function and derivatives
    R          = we*Ce*Be/wb;
    dRdxi(:,1) = we*Ce*(dBedxi(:,1)/wb-dwbdxi(1)*Be/(wb*wb));
    dRdxi(:,2) = we*Ce*(dBedxi(:,2)/wb-dwbdxi(2)*Be/(wb*wb));
    % Jacobian matrix
    dxdxi = pts'*dRdxi;
    %dxidx = inv(dxdxi);
    %dRdx  = dRdxi*dxidx;
    detJ  = det(dxdxi);
    mm    = R * R' * rho * detJ * wt;
    M(sctr,sctr) = M(sctr,sctr) + mm;
  end
end

for p=1:pCount
  body1.mass(p) = sum(M(p,:));
end

body2.mass = body1.mass;

bodies    = cell(1,1);
bodies{1} = body1;
bodies{2} = body2;
bodyCount = length(bodies);

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

for ib=1:length(bodies)
  body      = bodies{ib};
  elems     = ones(length(body.coord),1);
  
  for ip=1:length(body.coord)
    x = body.coord(ip,1); y = body.coord(ip,2);
    e = floor((x-xmin)/deltax) + 1 + numx2*floor((y-ymin)/deltay);
    elems(ip) = e;
  end
  
  bodies{ib}.elements = unique(elems);
  bodies{ib}.nodes    = unique(element(bodies{ib}.elements,:));
  mpoints = cell(elemCount,1);
  for ie=1:elemCount
    id  = find(elems==ie);
    mpoints{ie}=id;
  end
  
  bodies{ib}.mpoints  = mpoints;
end

%%
figure(1)
hold on
plot_mesh(node,element,'Q4','k-',1.);
nrbkntplot(circle1)
nrbctrlplot(circle1)
nrbkntplot(circle2)
nrbctrlplot(circle2)
view([0 90])
title('Initial configuration','FontSize',20)
set(gca,'FontSize',14)

%% grid data (mass, momentum, forces)

nmass      = zeros(nodeCount,1);  % nodal mass vector
nmomentum  = zeros(nodeCount,2);  % nodal momentum vector at the end of time step
nmomentum0 = zeros(nodeCount,2);  % nodal momentum vector at beginning of time step
nmomentumS = zeros(nodeCount,2);  % nodal momentum vector (mapped back in MUSL)
niforce    = zeros(nodeCount,2);  % nodal internal force vector
neforce    = zeros(nodeCount,2);  % nodal external force vector

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-14; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.001;
time  = 3.2;
t     = 0.;

interval = 100;
nsteps   = floor(time/dtime);

pos   = cell(nsteps/interval,2);

ta = [];           % time
ka = [];           % kinetic energy 
sa = [];           % strain energy

istep = 1;
ii = 1;
while ( t < time )
  disp(['time step ',num2str(t)])
  
  nmass(:)      = 0;
  nmomentum0(:) = 0;
  nmomentumS(:) = 0;
  niforce(:)    = 0;
  
  % loop over FEM bodies
  for ib=1:bodyCount
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints;
    for ie=1:length(elems)         % loop over computational cells or elements
      e     = elems(ie);
      esctr = element(e,:);      % element connectivity
      enode = node(esctr,:);     % element node coords
      mpts  = mpoints{e};        % particles inside element e
      for p=1:length(mpts)       % loop over particles
        pid    = mpts(p);
        xxp    = body.coord(pid,:);
        Mp     = body.mass(pid);
        vp     = body.velo(pid,:);
        forces = body.fp(pid,:);
        
        pt(1)= (2*xxp(1)-(enode(1,1)+enode(2,1)))/deltax;
        pt(2)= (2*xxp(2)-(enode(2,2)+enode(3,2)))/deltay;
        [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
        J0       = enode'*dNdxi;             % element Jacobian matrix
        %invJ0    = inv(J0);
        dNdx     = dNdxi/J0;
        % loop over nodes of current element "ie"
        for i=1:length(esctr)
          id               = esctr(i);
          Ni               = N(i);
          nmass(id)        = nmass(id)         + Ni*Mp;
          nmomentum0(id,:) = nmomentum0(id,:)  + Ni*Mp*vp;
          niforce(id,:)    = niforce(id,:)     - Ni*forces;
        end
      end
    end
  end
  % update nodal momenta
  nmomentum                 = nmomentum0 + niforce*dtime;
  %nmomentum(bGrid.tNodes,:) = 0;
  %nmomentum(bGrid.bNodes,:) = 0;
  
  %% update particle velocity and positions for all particles
  k = 0;
  for ib=1:bodyCount
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints;
    % loop over computational cells or elements
    for ie=1:length(elems)
      e     = elems(ie);
      esctr = element(e,:);      % element connectivity
      enode = node(esctr,:);     % element node coords
      mpts  = mpoints{e};        % particles inside element e
      % loop over particles
      for p=1:length(mpts)
        pid  = mpts(p);
        xp   = body.coord(pid,:);
        xp0  = body.coord0(pid,:);
        Mp   = body.mass(pid);
        vp   = body.velo(pid,:);
        % retrieve the grid functions/grads
        %N    = gBasis(ib,pid,:);
        for i=1:length(esctr)
          id      = esctr(i);
          x       = xp0 - node(id,:);
          [N,~]   = getMPM2D(x,deltax,deltay);
          if nmass(id) > 0
            massInv  = (1/nmass(id))*N;
            vp       = vp +  (nmomentum(id,:)-nmomentum0(id,:))*massInv;
            xp       = xp + dtime * nmomentum(id,:)*massInv;
          end
        end
        bodies{ib}.velo(pid,:)  = vp;
        bodies{ib}.coord(pid,:) = xp;
        
        k = k + 0.5*(vp(1)^2+vp(2)^2)*Mp;
        
        % mapped back bGrid momenta (used to compute L,,epsilon and stress)
        for i=1:length(esctr)
          id      = esctr(i);
          x       = xp0 - node(id,:);
          [N,~]   = getMPM2D(x,deltax,deltay);
          nmomentumS(id,:)  = nmomentumS(id,:) +  Mp*vp*N;
        end
      end
    end
    % as control points moved, build updated NURBS (for visualization)
    nurbs = buildNURBS(bodies{ib}.nurbs.knots,bodies{ib}.coord,...
        bodies{ib}.mesh.weights,bodies{ib}.mesh.noPtsX,bodies{ib}.mesh.noPtsY);
    bodies{ib}.nurbs     = nurbs;  
    bodies{ib}.mesh.node = bodies{ib}.coord;  % don;t forget to update mesh.node
  end
  %nmomentumS(bGrid.tNodes,:) = 0;
  %nmomentumS(bGrid.bNodes,:) = 0;
  % update  particles forces (step1)
  for ib=1:bodyCount
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints;
    up        = body.up; up(:) = 0.;
    % loop over computational cells or elements
    for ie=1:length(elems)
      e     = elems(ie);
      esctr = element(e,:);      % element connectivity
      enode = node(esctr,:);     % element node coords
      mpts  = mpoints{e};        % particles inside element e
      % loop over particles
      for p=1:length(mpts)
        pid  = mpts(p);
        xp0  = body.coord0(pid,:);
        % retrieve the grid functions/grads
        %N    = gBasis(ib,pid,:);
        for i=1:length(esctr)
          id = esctr(i);
           x      = xp0 - node(id,:);
          [N,~]   = getMPM2D(x,deltax,deltay);
          if nmass(id) > 0
            vI  = nmomentum(id,:)/nmass(id);   % nodal velocity
            up(pid,:) = up(pid,:) + dtime*vI*N; % particle increment displacement
          end
        end
      end
    end
    bodies{ib}.up = up;
  end
  % update FEM particles forces(step2)
  u = 0;
  for ib=1:bodyCount
    mesh = bodies{ib}.mesh;
    up   = bodies{ib}.up;
    f    = zeros(mesh.nodeCount*2,1);
    for e=1:mesh.elemCount
      sctr   = mesh.element(e,:);         %  element scatter vector
      nn     = length(sctr);
      sctrB(1,1:nn)      = sctr;
      sctrB(1,nn+1:2*nn) = sctr+mesh.nodeCount;
      B      = zeros(3,2*nn);
      Ce     = mesh.bezierExtractor(:,:,e); % element Bezier extraction operator
      we     = diag(mesh.weights(sctr));    % element weights
      pts    = mesh.node(sctr,:);           % element nodes' coords
      Wb     = Ce'*mesh.weights(sctr);      % element Bezier weights
      deltaU = up(sctr,:);                  % element displacements matrix nn x 2
      deltaU = deltaU(:);                   % make it a column vector 2nn x 1
      for gp=1:size(W,1)                    % loop over Gauss points
        pt      = Q(gp,:);
        wt      = W(gp);
        % Bernstein basis and derivatives at GP gp
        Be      = shapes(gp,:)';
        dBedxi  = reshape(derivs(gp,:,:),noBasis,2);
        % Bezier weight functions (denomenator of NURBS)
        wb        = dot(Be,Wb);
        dwbdxi(1) = dot(dBedxi(:,1),Wb);
        dwbdxi(2) = dot(dBedxi(:,2),Wb);
        % Shape function and derivatives
        R          = we*Ce*Be/wb;
        dRdxi(:,1) = we*Ce*(dBedxi(:,1)/wb-dwbdxi(1)*Be/(wb*wb));
        dRdxi(:,2) = we*Ce*(dBedxi(:,2)/wb-dwbdxi(2)*Be/(wb*wb));
        % Jacobian matrix
        dxdxi = pts'*dRdxi;
        %dxidx = inv(dxdxi);
        dRdx  = dRdxi/dxdxi;
        detJ  = det(dxdxi);
        % B matrix
        B(1,1:nn)      = dRdx(:,1)';
        B(2,nn+1:2*nn) = dRdx(:,2)';
        B(3,1:nn)      = dRdx(:,2)';
        B(3,nn+1:2*nn) = dRdx(:,1)';
        dEps           = B * deltaU;
        dSig           = C * dEps;
        bodies{ib}.strain(e,gp,:) = squeeze(bodies{ib}.strain(e,gp,:)) + dEps;
        bodies{ib}.stress(e,gp,:) = squeeze(bodies{ib}.stress(e,gp,:)) + dSig;
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        f(sctrB) = f(sctrB) + B' *squeeze(bodies{ib}.stress(e,gp,:)) * detJ * wt;
        
        u = u + 0.5*detJ*wt*dot(squeeze(bodies{ib}.stress(e,gp,:)),...
                                squeeze(bodies{ib}.strain(e,gp,:)));
      end
    end
    bodies{ib}.fp = reshape(f,mesh.nodeCount,2);
  end
  
  % update the element particle list
  for ib=1:bodyCount
    body      = bodies{ib};
    elems     = ones(length(body.volume),1);
    bodies{ib}.coord0 = bodies{ib}.coord;
    for ip=1:length(body.volume)
      x = body.coord(ip,1);
      y = body.coord(ip,2);
      m = body.mass(ip);
      e = floor((x-xmin)/deltax) + 1 + numx2*floor((y-ymin)/deltay);
      elems(ip) = e;
    end
    bodies{ib}.elements = unique(elems);
    bodies{ib}.nodes    = unique(element(bodies{ib}.elements,:));
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
  
  % VTK output
  
  if (  mod(istep-1,interval) == 0 )
    pos{ii,1} = bodies{1}.nurbs;
    pos{ii,2} = bodies{2}.nurbs; ii = ii + 1;
  end
  
  % advance to the next time step
  t = t + dtime;
  istep = istep + 1;
end

%%
figure(2)
hold on
plot_mesh(node,element,'Q4','k-',1.);
nrbctrlplot(bodies{1}.nurbs)
nrbctrlplot(bodies{2}.nurbs)
%nrbkntplot(circle1)
%nrbkntplot(circle2)
view([0 90])
title('Deformed configuration','FontSize',20)
set(gca,'FontSize',14)

%%
figure 
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(ta(1:end),sa(1:end),'r--','LineWidth',2);
plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy')
legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%axis([0 3 0 3])
grid on
box on
kinetic = 2*0.5*rho*pi*radius^2*(v0^2+v0^2);
%%
v = VideoWriter('iga-mpm.avi');
open(v);
for i=1:length(pos)
  figure
  hold on
  plot_mesh(node,element,'Q4','k-',1.);
  n1 = pos{i,1}; n2 = pos{i,2};
  nrbkntplot(n1)
  nrbkntplot(n2)
  view([0 90])
  axis off
  title( sprintf('%i',num2str(90*dtime*i) ));
  frame(i) = getframe;
  writeVideo(v,frame(i));
end
close all;
movie(frame,10);
close(v);
