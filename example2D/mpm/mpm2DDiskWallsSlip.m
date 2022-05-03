% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
%
%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 25 August 2015.

%%

addpath ../../grid/
addpath ../../basis/
addpath ../../particleGen/
addpath ../../constitutiveModels/
addpath ../../util/
addpath ../../postProcessing/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%

E   = 0.073e9;     % Young's modulus
nu  = 0.4;         % Poisson ratio
rho = 1.01e3;        % density
mu    = E/2/(1+nu);% shear modulus

v0    = 20;      % initial particle velocity [m/s]

interval      = 20;% time interval for saving vtp files.
vtkFileName   = 'mpm2DDiskWallSlip';
vtkFileName1  = '../../results/mpm/mpmDiskWallsGrid';

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

contact = 1;
plotNormal=0;

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid
lx = 12e-2;
ly = 12e-2;
ghostCell=0;
numx2 = 40;      % number of elements along X direction
numy2 = 40;      % number of elements along Y direction
[mesh]= buildGrid2D(lx,ly,numx2,numy2, ghostCell);
%mesh.node = mesh.node - [ones(mesh.nodeCount,1)*mesh.deltax ones(mesh.nodeCount,1)*mesh.deltay];

element   = mesh.element;
node      = mesh.node;
elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
deltax    = mesh.deltax;
deltay    = mesh.deltay;


% find boundary nodes

eps=1e-12;
left   = mesh.lNodes;
right  = mesh.rNodes;

%%   particle distribution

circle1.center = (1/100)*[6 6];
circle2.center = circle1.center;
circle1.radius = 3e-2; % m
circle2.radius = 4e-2; % m

ppc = [2 2];

[res] = generateMPDiffCircles(circle1,circle2,ppc,mesh);
pCount = size(res.position,1);
body1.volume = res.volume;
body1.volume0 = res.volume;
body1.mass   = res.volume*rho;
body1.coord  = res.position;
body1.deform = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress = zeros(pCount,3);                % stress
body1.strain = zeros(pCount,3);                % strain
body1.velo   = zeros(pCount,2);               % velocity
body1.velo(:,1) = v0;

% rigid particles

ppc = [3 3];

quad = [lx-mesh.deltax 0; lx 0; lx ly; lx-mesh.deltax ly];

[res] = generateMPForQuad(quad,ppc,mesh);
pCount       = size(res.position,1);
body2.volume = res.volume;
body2.volume0 = res.volume;
body2.mass   = res.volume*1e8;
body2.coord  = res.position;
body2.deform = repmat([1 0 0 1],pCount,1);     % gradient deformation
body2.stress = zeros(pCount,3);                % stress
body2.strain = zeros(pCount,3);                % strain
body2.velo   = zeros(pCount,2);               % velocity

quad = [0 0; mesh.deltax 0; mesh.deltax ly; 0 ly];

[res] = generateMPForQuad(quad,ppc,mesh);
pCount       = size(res.position,1);
body3.volume = res.volume;
body3.volume0 = res.volume;
body3.mass   = res.volume*rho;
body3.coord  = res.position;
body3.deform = repmat([1 0 0 1],pCount,1);     % gradient deformation
body3.stress = zeros(pCount,3);                % stress
body3.strain = zeros(pCount,3);                % strain
body3.velo   = zeros(pCount,2);               % velocity

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

% store bodies in this order: deformable bodies first, followed by rigid
% bodies
bodies(1)  = body1;
bodies(2)  = body2;
bodies(3)  = body3;

dBodyCount = 1;
bodyCount  = length(bodies);

for ib=1:length(bodies)
  body      = bodies(ib);
  elems     = ones(length(body.volume),1);
  
  for ip=1:length(body.volume)
    x = body.coord(ip,1); y = body.coord(ip,2);
    e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
    elems(ip) = e;
  end
  
  bodies(ib).elements = unique(elems);
  bodies(ib).nodes    = unique(element(bodies(ib).elements,:));
  
  mpoints = cell(elemCount,1);
  for ie=1:elemCount
    id  = find(elems==ie);
    mpoints{ie}=id;
  end
  
  bodies(ib).mpoints  = mpoints;
end


%% plot mesh, particles

figure
hold on
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(bodies{1}.elements,:),'Q4','cy-',2.1);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','r-',2.1);
plot(bodies(1).coord(:,1),bodies(1).coord(:,2),'k.','markersize',10);
plot(bodies(2).coord(:,1),bodies(2).coord(:,2),'b.','markersize',10);
plot(bodies(3).coord(:,1),bodies(3).coord(:,2),'b.','markersize',10);
axis off

%% check grid normals in a multiple body case

figure(100)
hold on
plot_mesh(node,element,'Q4','k-',1.);
plot(body1.coord(:,1),body1.coord(:,2),'k.','markersize',10);
plot(body2.coord(:,1),body2.coord(:,2),'b.','markersize',10);
plot(body3.coord(:,1),body3.coord(:,2),'b.','markersize',10);


for ib=1:bodyCount
  body      = bodies(ib);
  nodes     = body.nodes;
  [cellDensity,normals] = computeGridNormal(mesh,body);
  
  for i=1:length(nodes)
    nid = nodes(i);
    xI  = node(nid,:);
    nI  = normals(nid,:);
    nI = nI / norm(nI);
    le = 0.06;
    plot([xI(1) xI(1)+le*nI(1)],[xI(2) xI(2)+le*nI(2)],'r-','LineWidth',2);
  end
end
axis off


%% node quantities

nmass      = zeros(nodeCount,bodyCount+1);      % nodal mass vector
nmomentum  = zeros(nodeCount,2*(bodyCount+1));  % nodal momentum vector at the end
nmomentum0 = zeros(nodeCount,2*(bodyCount+1));  % nodal momentum vector at the begin
niforce    = zeros(nodeCount,2*(bodyCount+1));

%%
disp([num2str(toc),'   SOLVING '])

tol   = 1e-12; % mass tolerance

c     = sqrt(E/rho);
dtime = 1e-6;
time  = 1e-2; %time=9e-06;
t     = 0;

nsteps = floor(time/dtime);

istep = 1;

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

while ( t < time )
  disp(['time step ',num2str(t)])
  
  niforce(:)    = 0;
  nmass(:)      = 0;
  nmomentum0(:) = 0;
  nmomentum(:)  = 0;
  % loop over deformable bodies (update nodal momenta without contact)
  for ib = 1 : dBodyCount
    body      = bodies(ib);
    elems     = body.elements;
    mpoints   = body.mpoints;
    bodyId    = 2*ib-1:2*ib;       % index of body 'ib' in the nodal quantities
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
        stress =  bodies(ib).stress(pid,:);
        
        pt(1)= (2*xp(1)-(enode(1,1)+enode(2,1)))/deltax;
        pt(2)= (2*xp(2)-(enode(2,2)+enode(3,2)))/deltay;
        [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
        J0       = enode'*dNdxi;             % element Jacobian matrix
        invJ0    = inv(J0);
        dNdx     = dNdxi*invJ0;
        % loop over nodes of current element "ie"
        for i=1:length(esctr)
          id    = esctr(i);
          dNIdx = dNdx(i,1);
          dNIdy = dNdx(i,2);
          nmass(id,ib)          = nmass(id,ib)           + N(i)*Mp;
          nmomentum0(id,bodyId) = nmomentum0(id,bodyId)  + N(i)*Mp*vp;
          niforce(id,bodyId(1)) = niforce(id,bodyId(1)) - Vp*(stress(1)*dNIdx + stress(3)*dNIdy);
          niforce(id,bodyId(2)) = niforce(id,bodyId(2)) - Vp*(stress(3)*dNIdx + stress(2)*dNIdy);
        end
      end     % end loop over particles
    end       % end loop over cells of body 'ib'
    
    % update momenta
    nmomentum(:,bodyId) = nmomentum0(:,bodyId) + niforce(:,bodyId)*dtime;
    
    % compute system nodal quantities
    sysId                = 2*(bodyCount+1)-1:2*(bodyCount+1);
    nmomentum(:,sysId)   = nmomentum(:,sysId)   + nmomentum(:,bodyId);
    nmass(:,bodyCount+1) = nmass(:,bodyCount+1) + nmass(:,ib);
  end           % end loop over bodies
  
  % find contact nodes, actually common nodes between two bodies
  % loop over rigid bodies
  for irb = dBodyCount+1:bodyCount
    body      = bodies(irb);
    nodes     = body.nodes;
    % the following assumed contact between body1 and one of the rigid
    % bodies!!! @todo: generalise it
    contactNodes = intersect(bodies(1).nodes,nodes);
    if ~isempty(contactNodes)
      disp('contact is happening')
      [cellDensity,normals] = computeGridNormal(mesh,body);
      bodies(1).normals     = normals;
      
      % correct contact node velocities for deformable bodies
      for ib=1:dBodyCount
        for in=1:length(contactNodes)
          id       = contactNodes(in);
          velo1    = nmomentum(id,2*ib-1:2*ib);
          
          velocm   = [0 0]; % velocity of stationary rigid wall
          nI       = -bodies(ib).normals(id,:);
          nI       = nI / norm(nI);
          alpha    = dot(velo1 - velocm, nI); % shoud be: velo1-nmass*velocm
          if ( alpha > 0 )
            nmomentum(id,2*ib-1:2*ib) = velo1 - alpha*nI;
            disp('approaching')
          else
            disp('separating')
          end
        end
      end
    end
  end
  
  % update particle velocity  of deformable bodies
  
  for ib=1:dBodyCount
    body      = bodies(ib);
    elems     = body.elements;
    mpoints   = body.mpoints;
    bodyId    = 2*ib-1:2*ib;
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
        
        pt(1)= (2*xp(1)-(enode(1,1)+enode(2,1)))/deltax;
        pt(2)= (2*xp(2)-(enode(2,2)+enode(3,2)))/deltay;
        
        [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
        
        for i=1:length(esctr)
          id = esctr(i);
          aI = (nmomentum(id,bodyId)-nmomentum0(id,bodyId))/nmass(id,ib);
          vp = vp  + N(i)*aI;
        end
        bodies(ib).velo(pid,:) = vp;
      end
    end
  end
  
  % compute velocity to use in strain increments
  % mapped back from updated particle velocity
  
  for ib=1:dBodyCount
    body      = bodies(ib);
    elems     = body.elements;
    mpoints   = body.mpoints;
    bodyId    = 2*ib-1:2*ib;
    nmomentum(:,bodyId) = 0;
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
        
        pt(1)= (2*xp(1)-(enode(1,1)+enode(2,1)))/deltax;
        pt(2)= (2*xp(2)-(enode(2,2)+enode(3,2)))/deltay;
        
        [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
        
        for i=1:length(esctr)
          id = esctr(i);
          nmomentum(id,bodyId) = nmomentum(id,bodyId) + N(i)*Mp*vp;
        end
      end
    end
    
    % apply BCs on the mapped back grid momenta
    if ~isempty(contactNodes)
      % correct contact node velocities for deformable bodies
      
      for in=1:length(contactNodes)
        id       = contactNodes(in);
        velo1    = nmomentum(id,2*ib-1:2*ib);
        
        velocm   = [0 0]; % velocity of stationary rigid wall
        nI       = -bodies(ib).normals(id,:);
        nI       = nI / norm(nI);
        alpha    = dot(velo1 - velocm, nI); % shoud be: velo1-nmass*velocm
        if ( alpha > 0 )
          nmomentum(id,2*ib-1:2*ib) = velo1 - alpha*nI;
          %disp('approaching')
        else
          %disp('separating')
        end
      end
    end
    
  end
  
  % update particle position, stress, ....
  
  k = 0; u = 0;
  for ib=1:dBodyCount
    body      = bodies(ib);
    elems     = body.elements;
    mpoints   = body.mpoints;
    bodyId    = 2*ib-1:2*ib;
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
        Vp   = body.volume(pid);
        
        pt(1)= (2*xp(1)-(enode(1,1)+enode(2,1)))/deltax;
        pt(2)= (2*xp(2)-(enode(2,2)+enode(3,2)))/deltay;
        
        [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
        J0       = enode'*dNdxi;             % element Jacobian matrix
        invJ0    = inv(J0);
        dNdx     = dNdxi*invJ0;
        
        Lp   = zeros(2,2);
        for i=1:length(esctr)
          id = esctr(i); vI = [0 0];
          if nmass(id,ib) > 0
            vI = nmomentum(id,bodyId)/nmass(id,ib);
          end
          xp  = xp  + dtime * N(i)*vI;
          Lp  = Lp  + vI'*dNdx(i,:);
        end
        
        bodies(ib).coord(pid,:)= xp;
        
        % update stress last
        
        F       = ([1 0;0 1] + Lp*dtime)*reshape(bodies(ib).deform(pid,:),2,2);
        bodies(ib).deform(pid,:) = reshape(F,1,4);
        bodies(ib).volume(pid  ) = det(F)*bodies(ib).volume0(pid);
        dEps    = dtime * 0.5 * (Lp+Lp');
        dsigma  = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
        bodies(ib).stress(pid,:)  = bodies(ib).stress(pid,:) + dsigma';
        bodies(ib).strain(pid,:)  = bodies(ib).strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
        
        % compute strain, kinetic energies
        vp   = bodies(ib).velo(pid,:);
        k = k + 0.5*(vp(1)^2+vp(2)^2)*Mp;
        u = u + 0.5*Vp*bodies(ib).stress(pid,:)*bodies(ib).strain(pid,:)';
      end
    end
  end
  
  % update the element particle list
  
  for ib=1:dBodyCount
    body      = bodies(ib);
    elems     = ones(length(body.volume),1);
    
    for ip=1:length(body.volume)
      x = body.coord(ip,1); y = body.coord(ip,2);
      e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
      elems(ip) = e;
    end
    
    bodies(ib).elements = unique(elems);
    bodies(ib).nodes    = unique(element(bodies(ib).elements,:));
    mpoints = cell(elemCount,1);
    for ie=1:elemCount
      id  = find(elems==ie);
      mpoints{ie}=id;
    end
    
    bodies(ib).mpoints  = mpoints;
  end
  
  % store time,velocty for plotting
  
  ta = [ta;t];
  ka = [ka;k];
  sa = [sa;u];
  
  % VTK output
  
  if (  mod(istep-1,interval) == 0 )
    xp = [bodies(1).coord;bodies(2).coord;bodies(3).coord];
    s  = [bodies(1).stress;bodies(2).stress;bodies(3).stress];
    stress = [s sum(s,2)/3];
    data.stress = stress;
    vtkFile = sprintf('../../results/mpm/diskWall/%s%d',vtkFileName,istep-1);
    VTKParticles(xp,vtkFile,data);
  end
  
  
  % advance to the next time step
  
  t = t + dtime;
  istep = istep + 1;
end

%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

%%


% pvdFile = fopen(strcat('../results/',vtkFileName,'.pvd'), 'wt');
%
% fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
% fprintf(pvdFile,'<Collection>\n');
%
% for i = 1:nsteps
%      if (  mod(i,interval) == 0 )
%     vtuFile = sprintf('%s%d%s',vtkFileName,i,'.vtp');
%     fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
%      end
% end
%
% fprintf(pvdFile,'</Collection>\n');
% fprintf(pvdFile,'</VTKFile>\n');
%
% fclose(pvdFile);


Ux= zeros(size(node,1),1);
Uy= zeros(size(node,1),1);
sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

VTKPostProcess(node,element,2,'Quad4',vtkFileName1,...
    [sigmaXX sigmaYY sigmaXY],[Ux Uy]);
  
%ka = ka*1000;
%sa = sa*1000;

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(ta(1:end),sa(1:end),'r--','LineWidth',2);
plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
% xlabel('Time')
% ylabel('Energy (x1E-3)')
% legend('kinetic','strain','total')
% %set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
% %axis([0 3 0 3])


% savefile = 'mpm-explicit-2disks-dt001.mat';
% save(savefile,'ta','ka','sa');

disp([num2str(toc),'   DONE '])

