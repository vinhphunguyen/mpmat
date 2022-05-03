% This file implements the  coupled FEM-Material Point Method
%
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
%
% Disk-wire problem by York IJNME 1999.
%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 14 September 2015.

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
clear 
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

mkdir('../../results/mpm/diskWire');
vtkFileName   = '../../results/mpm/diskWire/mpmDiskWire';
vtkFileName1  = '../../results/mpm/diskWire/grid';

contact       = 1; %either 1 or 0 (free contact in MPM)

%% Material properties

rhoD     = 1;
rhoW     = 0.5;
youngD   = 1e4;
youngW   = 1e4;
poissonD = 0.3;
poissonW = 0.0;

L0       = 10-2e-3;   %% TODO: this is to ensure particle not exactly on the grid edge
v0       = 3.0; 
A        = 1;

stressState = 'PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C           = elasticityMatrix(youngD,poissonD,stressState);

tic;

%% Computational grid (prefer square cells)

l = 12;
w = 10;      % 10: length of the wire

noY0      = 30;         % number of elements along Y direction
noX0      = noY0*(l/w)-1; % number of elements along X direction
ghostCell = 0;

[bGrid]    = buildGrid2D(l,w,noX0,noY0, ghostCell);

node      = bGrid.node;
element   = bGrid.element;
deltax    = bGrid.deltax;
deltay    = bGrid.deltay;
elemCount = bGrid.elemCount;
nodeCount = bGrid.nodeCount;
numx2     = bGrid.numx;
numy2     = bGrid.numy;
Vc        = deltax * deltay;

%% generate material points
% disk 1 and disk 2
ppc           = [2 2];
disk1.center  = [3.5 6]; % left disk
disk2.center  = [8.5 4]; % right disk
disk1.radius  = 1.5;
disk2.radius  = 1.5;

[res1]         = generateMPForCircle(disk1,ppc,bGrid);
[res2]         = generateMPForCircle(disk2,ppc,bGrid);

pCount         = size(res1.position,1);
body1.volume   = res1.volume;
body1.volume0  = res1.volume;
body1.mass     = res1.volume*rhoD;
body1.coord    = res1.position;
body1.coord0   = res1.position;
body1.deform   = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress   = zeros(pCount,3);                % stress
body1.strain   = zeros(pCount,3);                % strain
body1.velo     = zeros(pCount,2);                % velocity
body1.velo(:,1) = v0;
body1.color    = ones(pCount,1);

pCount         = size(res2.position,1);
body2.volume   = res2.volume;
body2.volume0  = res2.volume;
body2.mass     = res2.volume*rhoD;
body2.coord    = res2.position;
body2.coord0   = res2.position;
body2.deform   = repmat([1 0 0 1],pCount,1);     % gradient deformation
body2.stress   = zeros(pCount,3);                % stress
body2.strain   = zeros(pCount,3);                % strain
body2.velo     = zeros(pCount,2);                % velocity
body2.velo(:,1) = -v0;
body2.color    = 2*ones(pCount,1);

% wire
pCount     = 30;
ghostCell  = 0;
[wGrid]    = buildGrid1D(L0,pCount-1, ghostCell);

Mp  = zeros(pCount,1);                % mass
xp  = zeros(pCount,2);                % position
vp  = zeros(pCount,2);                % velocity
fp  = zeros(pCount,2);                % particle forces (due to FEM-MPM)
up  = zeros(pCount,2);                % particle displacement (due to FEM-MPM)
le  = zeros(pCount-1,1);              % element length (due to FEM-MPM)
estress = zeros(pCount-1,1);           % element stress (due to FEM-MPM)
estrain = zeros(pCount-1,1);           % element stress (due to FEM-MPM)
xp(:,1) = l/2;                        % wire at the middle of the domain
xp(:,2) = wGrid.node + 1e-3;
Mp(:)   = rhoW*L0/pCount;
Mp(1)   = Mp(2)/2;
Mp(end) = Mp(2)/2;

wire.coord  = xp;
wire.coord0 = xp;
wire.mass   = Mp;
wire.velo   = vp;
wire.fp     = fp;
wire.up     = up;
wire.volume = Mp;
wire.color  = 3*ones(pCount,1);

% put all MPM bodies in one variable
bodies    = cell(3,1);
bodies{1} = body1;
bodies{2} = body2;
bodies{3} = wire;
bodyCount = length(bodies);
mpmBodyCount = bodyCount-1;
%%
figure(1)
set(gca,'FontSize',14)
hold on
plot_mesh(node,element,'Q4','k-',1.);
plot(body1.coord(:,1),body1.coord(:,2),'k.','markersize',22);
plot(body2.coord(:,1),body2.coord(:,2),'r.','markersize',22);
%plot(xp(:,1),xp(:,2),'b.','markersize',20,'MarkerFaceColor','b','LineWidth',2);
plot_mesh(xp,wGrid.element,'L2','r-*',2)
title('Initial configuration','FontSize',20)
set(gca,'FontSize',14)

%clear body1 body2  

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

for ib=1:length(bodies)
  body      = bodies{ib};
  elems     = ones(length(body.coord),1);
  
  for ip=1:length(body.coord)
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

%% grid data (mass, momentum, forces)

nmass      = zeros(nodeCount,1);  % nodal mass vector
nmomentum  = zeros(nodeCount,2);  % nodal momentum vector
nmomentum0 = zeros(nodeCount,2);  % nodal momentum vector
nmomentumS = zeros(nodeCount,2);  % nodal momentum vector
niforce    = zeros(nodeCount,2);  % nodal internal force vector
neforce    = zeros(nodeCount,2);  % nodal external force vector

% store grid basis/derivatives at particles of bodies
% only works if all bodies have the same number of particles!!!
% gBasis     = zeros(bodyCount,length(body1.volume),4);
% gGradX     = zeros(bodyCount,length(body1.volume),4);
% gGradY     = zeros(bodyCount,length(body1.volume),4);

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-14; % mass tolerance

c     = sqrt(youngW/rhoW);
dtime = 0.3*bGrid.deltax/c;
time  = 1.3;
t     = 0.;

interval = 100;
nsteps   = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

% kinetic and potential energies for disks and wire
ta  = [];
ke1 = [];
pe1 = [];
ke2 = [];
pe2 = [];
ke3 = [];
pe3 = [];

k1 = 0; p1 = 0; k2 = 0; p2 = 0; k3 = 0; p3 = 0;

while ( t < time )
  disp(['time step ',num2str(t)])
  
  nmass(:)      = 0;
  nmomentum0(:) = 0;
  nmomentumS(:) = 0;
  niforce(:)    = 0;

   
  % store time,velocty for plotting  
  ta  = [ta;t]; 
  ke1  = [ke1;k1];  ke2  = [ke2;k2]; ke3  = [ke3;k3];  
  pe1  = [pe1;p1];  pe2  = [pe2;p2]; pe3  = [pe3;p3]; 
  
  % loop over MPM bodies
  for ib=1:mpmBodyCount
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
        xxp     = body.coord(pid,:);
        Mp     = body.mass(pid);
        vp     = body.velo(pid,:);
        Vp     = body.volume(pid);
        stress =  bodies{ib}.stress(pid,:);
        
        pt(1)= (2*xxp(1)-(enode(1,1)+enode(2,1)))/deltax;
        pt(2)= (2*xxp(2)-(enode(2,2)+enode(3,2)))/deltay;
        [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
        J0       = enode'*dNdxi;             % element Jacobian matrix
        %invJ0    = inv(J0);
        dNdx     = dNdxi/J0; %dNdxi*invJ0;
        % store grid basis functions and gradients
%         gBasis(ib,pid,:) = N;
%         gGradX(ib,pid,:) = dNdx(:,1);
%         gGradY(ib,pid,:) = dNdx(:,2);
        % loop over nodes of current element "ie"
        for i=1:length(esctr)
          id    = esctr(i);
          dNIdx = dNdx(i,1);
          dNIdy = dNdx(i,2);
          nmass(id)       = nmass(id)        + N(i)*Mp;
          nmomentum0(id,:)= nmomentum0(id,:) + N(i)*Mp*vp;
          niforce(id,1)   = niforce(id,1) - Vp*(stress(1)*dNIdx + stress(3)*dNIdy);
          niforce(id,2)   = niforce(id,2) - Vp*(stress(3)*dNIdx + stress(2)*dNIdy);
        end
      end
    end
  end
  % loop over FEM bodies
  for ib=mpmBodyCount+1:bodyCount
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
        % store grid basis functions and gradients
        gBasis(ib,pid,:) = N;
        gGradX(ib,pid,:) = dNdx(:,1);
        gGradY(ib,pid,:) = dNdx(:,2);
        % loop over nodes of current element "ie"
        for i=1:length(esctr)
          id               = esctr(i);
          Ni               = N(i);
          nmass(id)        = nmass(id)         + Ni*Mp;
          nmomentum0(id,:) = nmomentum0(id,:)  + Ni*Mp*vp;
          niforce(id,:)    = niforce(id,:)     + Ni*forces;
        end
      end
    end
  end
  % update nodal momenta
  nmomentum                 = nmomentum0 + niforce*dtime;
  nmomentum(bGrid.tNodes,:) = 0;
  nmomentum(bGrid.bNodes,:) = 0;
  
  %% update particle velocity and positions for all particles (MPM and FEM)
  for ib=1:bodyCount
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints;
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
        xp0  = body.coord0(pid,:);
        Mp   = body.mass(pid);
        vp   = body.velo(pid,:);        
        % retrieve the grid functions/grads
        %N    = gBasis(ib,pid,:);
        for i=1:length(esctr)
          id      = esctr(i);
          x       = xp0 - node(id,:);
          [N,~]   = getMPM2D(x,deltax,deltay);
          if nmass(id) > tol
            massInv  = (1/nmass(id))*N;
            vp       = vp +  (nmomentum(id,:)-nmomentum0(id,:))*massInv;
            xp       = xp + dtime * nmomentum(id,:)*massInv;
          end
        end        
        bodies{ib}.velo(pid,:)  = vp;
        bodies{ib}.coord(pid,:) = xp;
        % mapped back bGrid momenta (used to compute L,,epsilon and stress)
        for i=1:length(esctr)
          id      = esctr(i);          
          x       = xp0 - node(id,:);
          [N,~]   = getMPM2D(x,deltax,deltay);
          nmomentumS(id,:)  = nmomentumS(id,:) +  Mp*vp*N;
        end
      end
    end
  end

  k1 = 0.5*dot( bodies{1}.mass, bodies{1}.velo(:,1).^2+bodies{1}.velo(:,2).^2);
  k2 = 0.5*dot( bodies{2}.mass, bodies{2}.velo(:,1).^2+bodies{2}.velo(:,2).^2);
  k3 = 0.5*dot( bodies{3}.mass, bodies{3}.velo(:,1).^2+bodies{3}.velo(:,2).^2);

  nmomentumS(bGrid.tNodes,:) = 0;
  nmomentumS(bGrid.bNodes,:) = 0;
  
  bodies{2}.coord(1,:)   = wire.coord(1,:);
  bodies{2}.coord(end,:) = wire.coord(end,:);
  bodies{2}.velo(1,:)    = [0 0];
  bodies{2}.velo(end,:)  = [0 0];
  
  % update MPM particle stress
  for ib=1:mpmBodyCount
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints; 
    % loop over computational cells or elements
    for ie=1:length(elems)
      e     = elems(ie);
      esctr = element(e,:);      % element connectivity
      enode = node(esctr,:);     % element node coords
      mpts  = mpoints{e};       % particles inside element e
      % loop over particles
      for p=1:length(mpts)
        pid  = mpts(p);             
        xp0  = body.coord0(pid,:);
        % retrieve the grid functions/grads        
        %dNdx = [squeeze(gGradX(ib,pid,:)) squeeze(gGradY(ib,pid,:))];        
        Lp   = zeros(2,2);
        for i=1:length(esctr)
          id        = esctr(i);
          x         = xp0 - node(id,:);
          [~,dNdx]  = getMPM2D(x,deltax,deltay);
          if nmass(id) > 0
            vI  = nmomentumS(id,:)/nmass(id);  % nodal velocity            
            Lp  = Lp + vI'*dNdx;
          end
        end                          
        F       = ([1 0;0 1] + Lp*dtime)*reshape(bodies{ib}.deform(pid,:),2,2);
        bodies{ib}.deform(pid,:) = reshape(F,1,4);
        bodies{ib}.volume(pid  ) = det(F)*bodies{ib}.volume0(pid);
        dEps    = dtime * 0.5 * (Lp+Lp');
        dsigma  = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
        bodies{ib}.stress(pid,:)  = bodies{ib}.stress(pid,:) + dsigma';
        bodies{ib}.strain(pid,:)  = bodies{ib}.strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];        
       
        %u    = u + 0.5*Vp*bodies{ib}.stress(pid,:)*bodies{ib}.strain(pid,:)';
      end
    end
    bodies{ib}.coord0 = bodies{ib}.coord;
  end

  p1 = 0.5*dot( bodies{1}.volume,...
                dot( bodies{1}.stress,bodies{1}.strain,2 ) );


  p2 = 0.5*dot( bodies{2}.volume,...
                dot( bodies{2}.stress,bodies{2}.strain,2 ) );

  % update FEM particles forces, step1,
  for ib=mpmBodyCount+1:bodyCount
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
        % retrieve the grid functions/grads        
        N    = gBasis(ib,pid,:);       
        for i=1:length(esctr)
          id = esctr(i);
          if nmass(id) > 0
            vI  = nmomentumS(id,:)/nmass(id);  % nodal velocity            
            up(pid,:) = up(pid,:) + dtime*vI*N(i); % particle increment displacement
          end
        end                                 
      end
    end
    bodies{ib}.up = up;
  end
  % update FEM particles forces, assume only one set of FEM particles
  fp = bodies{3}.fp; fp(:)=0; p3 = 0;
  for e=1:size(wGrid.element,1)
      esctr    = wGrid.element(e,:);
      nodes    = bodies{bodyCount}.coord(esctr,:);
      ll       = norm(nodes(1,:)-nodes(2,:));
      cosTheta = (nodes(2,1)-nodes(1,1))/ll;
      sinTheta = (nodes(2,2)-nodes(1,2))/ll;
      Q        = [cosTheta  sinTheta;-sinTheta cosTheta];
      Qinv     = [cosTheta -sinTheta;sinTheta cosTheta];
      incDisp  = bodies{bodyCount}.up(esctr,:);
      deltaU1  = Q*incDisp(1,:)';
      deltaU2  = Q*incDisp(2,:)';            
      incStrain = (1/ll)*(-deltaU1(1)+deltaU2(1));      
      %if incStrain < 0, incStrain=0.; end
      %if incStrain < 0, disp('negative inc. strain'); end
      estrain(e) = estrain(e) + incStrain;
      sig        = estress(e) + youngW*incStrain;
      %if strain(e) < 0, disp('NEGATIVE STRAIN');sig=0;end      
      %if sig < 0, sig = 0; end
      estress(e) = sig;
      le(e)      = ll;
      force      = sig*A;
      fp(esctr(1),:)  = fp(esctr(1),:) + (Qinv*[force;0])';
      fp(esctr(2),:)  = fp(esctr(2),:) - (Qinv*[force;0])';
      p3 = p3 + estrain(e)*sig*0.5*ll;
  end
  bodies{3}.fp = fp;
  bodies{3}.stress = zeros(length(wire.coord),3);
  
  % update the element particle list  
  for ib=1:bodyCount
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

  
  % VTK output
  
  if (  mod(istep-1,interval) == 0 )
    xp = [bodies{1}.coord;bodies{2}.coord;bodies{3}.coord];
    s  = [bodies{1}.stress;bodies{2}.stress;bodies{3}.stress];
    stress = [s sum(s,2)/3];
    data.stress = stress;
    %data.velo=[bodies{1}.velo;bodies{2}.velo];
    data.color=[bodies{1}.color;bodies{2}.color;bodies{3}.color];
    vtkFile = sprintf('%s%d',vtkFileName,istep-1);
    VTKParticles(xp,vtkFile,data);
  end
    
  % advance to the next time step  
  t = t + dtime;
  istep = istep + 1;
end

%%

Ux= zeros(size(node,1),1);
Uy= zeros(size(node,1),1);
sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

VTKPostProcess(node,element,2,'Quad4',vtkFileName1,...
    [sigmaXX sigmaYY sigmaXY],[Ux Uy]);

%%
figure(2)
hold on
plot_mesh(node,element,'Q4','k-',1.);
plot(bodies{1}.coord(:,1),bodies{1}.coord(:,2),'k.','markersize',10);
plot(body1.coord(:,1),body1.coord(:,2),'r.','markersize',10);
plot(bodies{2}.coord(:,1),bodies{2}.coord(:,2),'k.','markersize',10);
plot(body2.coord(:,1),body2.coord(:,2),'r.','markersize',10);
plot(bodies{3}.coord(:,1),bodies{3}.coord(:,2),'bs','markersize',10,'MarkerFaceColor','b');
plot(node(bodies{3}.nodes,1),node(bodies{3}.nodes,2),'b*');
title('Deformed configuration','FontSize',20)
set(gca,'FontSize',14)

%%

% 
figure
subplot(2,2,1)
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ke1(1:end),'black-','LineWidth',2.1);
%plot(ta(1:end),ke2(1:end),'red-','LineWidth',2.1);
plot(ta(1:end),pe1(1:end),'red-','LineWidth',2.1);
%plot(ta(1:end),pe2(1:end),'red--','LineWidth',2.1);
%plot(ta(1:end),ke3(1:end),'blue-','LineWidth',2.1);
%plot(ta(1:end),ka(1:end)+sa(1:end),'blue-','LineWidth',2.1);
xlabel('Time')
ylabel('Energies')
legend('kinetic disk1','potential disk1')
grid on
box on
set(gca,'XTick',[0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4])

subplot(2,2,2)
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ke2(1:end),'black-','LineWidth',2.1);
%plot(ta(1:end),ke2(1:end),'red-','LineWidth',2.1);
plot(ta(1:end),pe2(1:end),'red-','LineWidth',2.1);
%plot(ta(1:end),pe2(1:end),'red--','LineWidth',2.1);
%plot(ta(1:end),ke3(1:end),'blue-','LineWidth',2.1);
%plot(ta(1:end),ka(1:end)+sa(1:end),'blue-','LineWidth',2.1);
xlabel('Time')
ylabel('Energies')
legend('kinetic disk2','potential disk2')
grid on
box on
set(gca,'XTick',[0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4])

subplot(2,2,3)
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ke3(1:end),'black-','LineWidth',2.1);
%plot(ta(1:end),ke2(1:end),'red-','LineWidth',2.1);
plot(ta(1:end),pe3(1:end),'red-','LineWidth',2.1);
%plot(ta(1:end),pe2(1:end),'red--','LineWidth',2.1);
%plot(ta(1:end),ke3(1:end),'blue-','LineWidth',2.1);
%plot(ta(1:end),ka(1:end)+sa(1:end),'blue-','LineWidth',2.1);
xlabel('Time')
ylabel('Energies')
legend('kinetic wire','potential wire')
grid on
box on
set(gca,'XTick',[0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4])

subplot(2,2,4)
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ke1(1:end)+pe1(1:end),'black-','LineWidth',2.1);
plot(ta(1:end),ke2(1:end)+pe2(1:end),'blue-','LineWidth',2.1);
plot(ta(1:end),ke3(1:end)+pe3(1:end),'red-','LineWidth',2.1);
plot(ta(1:end),ke1+pe1+ke2+pe2+ke3+pe3,'cyan-','LineWidth',2.9);
xlabel('Time')
ylabel('\Energies')
legend('total disk1','total disk2', 'total wire', 'total')
grid on
box on
set(gca,'XTick',[0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4])

