% This file implements the  Material Point Method
%
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
%
%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 1 September 2015.

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

I  = [1 0;0 1];

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


vtkFileName   = '../../results/mpm/square2Circle/mpmSquare2Circle';
vtkFileName1  = '../../results/mpm/square2Circle/grid';

contact       = 1; %either 1 or 0 (free contact in MPM)

%% Material properties
% dynes, cm, s

rho    = 1;
K_f    = 1.5e5;   % bulk modulus
gamma  = 7;       %
lambda = 0.5;     % dynamic viscosity
sigma  = 2.4;     % surface tension

tic;

%% Computational grid (all length in cm)

l = 2;
w = 2;

noX0      = 40;        % number of elements along X direction
noY0      = 40;        % number of elements along Y direction
ghostCell = 0;

[grid]    = buildGrid2D(l,w,noX0,noY0, ghostCell);

node      = grid.node;
element   = grid.element;
deltax    = grid.deltax;
deltay    = grid.deltay;
elemCount = grid.elemCount;
nodeCount = grid.nodeCount;
numx2     = grid.numx;
numy2     = grid.numy;
Vc        = deltax * deltay;

%% generate material points

ppc           = [2 2];
square        = [0.5 0.5; 1.5 0.5; 1.5 1.5; 0.5 1.5];

[res]          = generateMPForQuad(square,ppc,grid);

pCount         = size(res.position,1);
body1.volume   = res.volume;
body1.volume0  = res.volume;
body1.mass     = res.volume*rho;
body1.coord    = res.position;
body1.deform   = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress   = zeros(pCount,3);                % stress
body1.strain   = zeros(pCount,3);                % strain
body1.velo     = zeros(pCount,2);                % velocity
body1.pressure = zeros(pCount,1);                % pressure
body1.density  = rho*ones(pCount,1);             % density (change in time)


% put all bodies in one variable
bodies    = cell(1,1);
bodies{1} = body1;
bodyCount = length(bodies);

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

for ib=1:length(bodies)
  body      = bodies{ib};
  elems     = ones(length(body.volume),1);
  
  for ip=1:length(body.volume)
    x = body.coord(ip,1); y = body.coord(ip,2);
    e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
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

% boundary nodes

%% node quantities

nmassS    = zeros(nodeCount,1);  % nodal mass vector of the system
nmomentumS= zeros(nodeCount,2);  % nodal momentum vector of the system
niforceS  = zeros(nodeCount,2);  % nodal internal force of the system
nmass     = zeros(nodeCount,1);  % nodal mass vector of each body
nmomentum = zeros(nodeCount,2);  % nodal momentum vector of each body
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector
nvelo     = zeros(nodeCount,2*(bodyCount+1));  % nodal velocities (body1,body2,center of mass)
nvelo0    = nvelo;
nacce     = zeros(nodeCount,2*bodyCount);

gBasis     = zeros(bodyCount,pCount,4);
gGradX     = zeros(bodyCount,pCount,4);
gGradY     = zeros(bodyCount,pCount,4);


%% check grid normals in a multiple body case

figure(100)
hold on
coords=bodies{1}.coord;
plot_mesh(node,element,'Q4','k-',1.);
plot(coords(:,1),coords(:,2),'k.','markersize',10);

tol = 1e-14;
for ib=1:bodyCount
  body      = bodies{ib};
  nodes     = body.nodes;
  [res]     = computeGridNormalCurvature(grid,body);
  bodies{ib}.normals    = res.gNormals;
  bodies{ib}.curvatures = res.gCurvatures;
  normals   = res.gNormals;
  density   = res.cDensities;
  curvature = res.cCurvatures;
  for i=1:length(nodes)
    nid = nodes(i);
    xI  = node(nid,:);
    nI  = normals(nid,:);
    if norm(nI) > tol, nI = nI / norm(nI); end
    le = 0.1;
    %plot([xI(1) xI(1)+le*nI(1)],[xI(2) xI(2)+le*nI(2)],'r-','LineWidth',2);
    %quiver(xI(1), xI(2), le*nI(1), le*nI(2),4,'LineWidth',2,...
    %     'Color','r','MaxHeadSize',10);
    quiver(xI(1), xI(2), neforce(nid,1)/10, neforce(nid,2)/10,4,'LineWidth',2,...
         'Color','r','MaxHeadSize',2);       
  end
end
axis off
%%
figure
hold on
plot_field(node,element,'Q4',sqrt(sum(normals.^2,2)));
%plot_field(node,element,'Q4',curvature);
colormap jet
colorbar
h=colorbar;
set(h,'fontsize',14);
axis off

figure
hold on
patch('Faces',element,'Vertices',node,'FaceVertexCData',curvature,'FaceColor','flat')
plot(coords(:,1),coords(:,2),'k.','markersize',10);
colormap jet
colorbar
axis equal
axis off

ta  = [];           % time

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-14; % mass tolerance

c     = sqrt(K_f/rho);
dtime = 0.2*grid.deltax/c;
time  = 0.08;
t     = 0.;

interval     = 10;
nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
  disp(['time step ',num2str(t)])
  
  nvelo(:)     = 0;
  nmassS(:)     = 0;
  nmomentumS(:) = 0;
  %% loop over bodies (update nodal momenta without contact)
  for ib=1:bodyCount
    %% reset grid data (body contribution)
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    neforce(:)   = 0;
    
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints;
    jumpC     = rho;
    for ie=1:length(elems)         % loop over computational cells or elements
      e     = elems(ie);
      esctr = element(e,:);      % element connectivity
      enode = node(esctr,:);     % element node coords
      mpts  = mpoints{e};        % particles inside element e
      for p=1:length(mpts)       % loop over particles
        pid    = mpts(p);
        xp     = body.coord(pid,:);
        Mp     = body.mass(pid);
        vp     = body.velo(pid,:);
        Vp     = body.volume(pid);
        stress =  bodies{ib}.stress(pid,:);
        
        pt(1)= (2*xp(1)-(enode(1,1)+enode(2,1)))/deltax;
        pt(2)= (2*xp(2)-(enode(2,2)+enode(3,2)))/deltay;
        [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
        J0       = enode'*dNdxi;             % element Jacobian matrix
        invJ0    = inv(J0);
        dNdx     = dNdxi*invJ0;
        % store grid basis functions and gradients
        gBasis(ib,pid,:) = N;
        gGradX(ib,pid,:) = dNdx(:,1);
        gGradY(ib,pid,:) = dNdx(:,2);
        
        % loop over nodes of current element "ie"
        for i=1:length(esctr)
          id    = esctr(i);
          dNIdx = dNdx(i,1);
          dNIdy = dNdx(i,2);
          nmass(id)       = nmass(id)       + N(i)*Mp;
          nmomentum(id,:) = nmomentum(id,:) + N(i)*Mp*vp;
          niforce(id,1)   = niforce(id,1) - Vp*(stress(1)*dNIdx + stress(3)*dNIdy);
          niforce(id,2)   = niforce(id,2) - Vp*(stress(3)*dNIdx + stress(2)*dNIdy);
          % this is surface tension force
          neforce(id,:)   = (sigma/jumpC)*Vc*body.curvatures(id)*body.normals(id,:);
        end
      end
    end
    
    activeNodes = bodies{ib}.nodes;

    if (contact)
      massInv = 1./nmass(activeNodes);
      smallMassIds = find(massInv > 1e10);
      if ~isempty(smallMassIds)
        disp('small mass!!!')
        massInv(smallMassIds) = 0;
      end
      % store old velocity v_I^t
      nvelo0(activeNodes,2*ib-1) = nmomentum(activeNodes,1).*massInv;
      nvelo0(activeNodes,2*ib)   = nmomentum(activeNodes,2).*massInv;
      % update body nodal momenta
      nmomentum(activeNodes,:) = nmomentum(activeNodes,:) + (niforce(activeNodes,:) + neforce(activeNodes,:))*dtime;
      % store uncorrected updated body velocity and acceleration v_I^{t+\Delta t}
      nvelo(activeNodes,2*ib-1) = nmomentum(activeNodes,1).*massInv;
      nvelo(activeNodes,2*ib)   = nmomentum(activeNodes,2).*massInv;
      
      nacce(activeNodes,2*ib-1) = (niforce(activeNodes,1) + neforce(activeNodes,1)).*massInv;
      nacce(activeNodes,2*ib  ) = (niforce(activeNodes,2) + neforce(activeNodes,2)).*massInv;
    end
    
    % store system momentum and mass
    nmomentumS(activeNodes,:) = nmomentumS(activeNodes,:) + nmomentum(activeNodes,:);
    nmassS    (activeNodes  ) = nmassS    (activeNodes  ) + nmass(activeNodes);
    niforceS  (activeNodes,:) = niforceS  (activeNodes,:) + niforce(activeNodes,:);
  end
  

  % no boundary conditions 
    
  % find contact nodes, actually common nodes between two bodies
%   if contact
%     contactNodes = intersect(bodies{1}.nodes,bodies{2}.nodes);
%     if ~isempty(contactNodes)
%       %disp('contact is happening')
%       for ib=1:bodyCount
%         body      = bodies{ib};
%         nodes     = body.nodes;
%         [cellDensity,normals] = computeGridNormal(grid,body);
%         bodies{ib}.normals = normals;
%       end
%       % pause
%     end
%     
% 
%     %% correct contact node velocities
%     
%     for ib=1:bodyCount
%       for in=1:length(contactNodes)
%         id       =  contactNodes(in);
%         velo1    = nvelo(id,2*ib-1:2*ib);
%         velocm   = [0 0];
%         if nmassS(id) > tol
%           velocm   = nmomentumS(id,:)/nmassS(id);
%         else
%           disp('small mass detected')
%         end
%         nI       = bodies{ib}.normals(id,:);
%         nI       = nI / norm(nI);
%         deltaVe  = velo1 - velocm;
%         D        = dot(deltaVe, nI);
%         C        = deltaVe(1)*nI(2) - deltaVe(2)*nI(1);
%         absC     = abs(C);
%         muPrime  = min(fric,absC/D);
%         if ( D >= 0 )
%           nvelo(id,2*ib-1:2*ib) = velo1 - D*( nI + (muPrime/absC)*[nI(2)*C -nI(1)*C] );
%           nacce(id,2*ib-1:2*ib) = (1/dtime)*( nvelo(id,2*ib-1:2*ib) - nvelo0(id,2*ib-1:2*ib) );
%           % disp('approaching')
%         else
%           % disp('separating')
%         end
%       end
%     end
%   end
  
  %% update particle velocity and position and stresse
  for ib=1:bodyCount
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints;
    if (contact)
      indices = 2*ib-1:2*ib;
    else
      indices = 2*bodyCount+1:2*bodyCount+2;
    end
    % loop over computational cells or elements
    for ie=1:length(elems)
      e     = elems(ie);
      esctr = element(e,:);      % element connectivity
      enode = node(esctr,:);     % element node coords
      mpts  = mpoints{e};       % particles inside element e
      % loop over particles
      for p=1:length(mpts)
        pid  = mpts(p);
        xp   = body.coord(pid,:);
        Mp   = body.mass(pid);
        vp   = body.velo(pid,:);
        Vp   = body.volume(pid);
        % retrieve the grid functions/grads
        N    = gBasis(ib,pid,:);
        dNdx = [squeeze(gGradX(ib,pid,:)) squeeze(gGradY(ib,pid,:))];
        
        Lp   = zeros(2,2);
        for i=1:length(esctr)
          id = esctr(i);
          vI = nvelo(id,indices);
          aI = nacce(id,indices);
          vp  = vp  + dtime * N(i)*aI;
          xp  = xp  + dtime * N(i)*vI;
          Lp  = Lp + vI'*dNdx(i,:);
        end
        
        bodies{ib}.velo(pid,:) = vp;
        bodies{ib}.coord(pid,:)= xp;
        
        % update stress last
        
        F       = ([1 0;0 1] + Lp*dtime)*reshape(bodies{ib}.deform(pid,:),2,2);
        bodies{ib}.deform(pid,:) = reshape(F,1,4);
        bodies{ib}.volume(pid  ) = det(F)*bodies{ib}.volume0(pid);
        
        dEps          = 0.5 * (Lp+Lp');
        traceEps      =  ( dEps(1,1) + dEps(2,2) );
        density       = bodies{ib}.density(pid)/(1 + dtime*traceEps);
        pressure      = K_f * ( (density/rho)^gamma - 1 );
        bodies{ib}.stress(pid,:)  = -pressure*[1;1;0] + ...
           2    * lambda * [dEps(1,1); dEps(2,2); 2*dEps(1,2)] - ...
          (2/3) * lambda * traceEps * [1;1;0] ;
        bodies{ib}.density(pid)  = density;
        bodies{ib}.pressure(pid) = pressure;
      end
    end
  end
  
  % update the element particle list
  
  for ib=1:length(bodies)
    body      = bodies{ib};
    elems     = ones(length(body.volume),1);
    
    for ip=1:length(body.volume)
      x = body.coord(ip,1);
      y = body.coord(ip,2);
      m = body.mass(ip);
      e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
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
  
  % compute normals and curvatures at grid nodes for next step
  [res]     = computeGridNormalCurvature(grid,body);
  bodies{ib}.normals    = res.gNormals;
  bodies{ib}.curvatures = res.gCurvatures;
  % store time,velocty for plotting
  
  %ta  = [ta;t];
  
  % VTK output
  
%   if (  mod(istep-1,interval) == 0 )
%     xp = [bodies{1}.coord;bodies{2}.coord];
%     s  = [bodies{1}.stress;bodies{2}.stress];
%     stress = [s sum(s,2)/3];
%     data.stress = stress;
%     data.velo=[bodies{1}.velo;bodies{2}.velo];
%     vtkFile = sprintf('%s%d',vtkFileName,istep-1);
%     VTKParticles(xp,vtkFile,data);
%   end
  
  
  % advance to the next time step
  
  t = t + dtime;
  istep = istep + 1;
end


%
% Ux= zeros(size(node,1),1);
% Uy= zeros(size(node,1),1);
% sigmaXX = zeros(size(node,1),1);
% sigmaYY = zeros(size(node,1),1);
% sigmaXY = zeros(size(node,1),1);
%
% VTKPostProcess(node,element,2,'Quad4',vtkFileName1,...
%     [sigmaXX sigmaYY sigmaXY],[Ux Uy]);

%% post processing

disp([num2str(toc),'   POST-PROCESSING '])



% noslip = load('rolling-sphere-noslip.mat');
% wislip = load('rolling-sphere-slip.mat');

% figure
% set(gca,'FontSize',14)
% hold on
% plot(tt,xx1,'r-','LineWidth',2.2);
% plot(ta,xcm,'black--','LineWidth',2);
% %plot(noslip.ta,noslip.xcm,'r--','LineWidth',2);
% xlabel('Time [s]')
% ylabel('center of mass position [m]')
% legend('stick-analytical','stick-MPM')
% %set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
% %axis([0 2 0 16])

disp([num2str(toc),'   DONE '])
