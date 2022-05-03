% This file implements the  Material Point Method
%
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
%
% Failure of rock under impact. Damage model by Kipps, used by Das in SPH
% paper.
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 23 September 2015.
% Unsuccesfull attempt.

%%

addpath ../../grid/
addpath ../../basis/
addpath ../../particleGen/
addpath ../../constitutiveModels/
addpath ../../util/
addpath ../../geoMesh/
addpath ../../externals/PolyMesher/
addpath ../../postProcessing/

%%
clc
clear
colordef white

I         = [1 1 0];
oneThird  = 1./3.;
twoThird  = 2/3.;

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

mkdir('../../results/mpm/diskRock/')
vtkFileName   = '../../results/mpm/diskRock/diskRock';
vtkFileName1  = '../../results/mpm/diskRock/grid';



%% Material properties

bulkModulus   = 12.2e9; % bulk modulus, Pa
shearModulus  = 2.67e9; % shear modulus, Pa
density       = 2300;   % kg/m3
k             = 6.35e46;% Weibull parameter
m             = 12.8;   % Weibull parameter
v0            = 100;    % velocity, m/s
Cg            = 0.4*sqrt((bulkModulus)/density);
alpha         = 8*pi*Cg^3*k/(m+1)/(m+2)/(m+3);
mAlpha        = (m+3)*alpha^oneThird;
mDiv3         = m/3;
minDamage0    = 0.001;
maxDamage0    = 0.01;

tic;

%% Computational grid (all length in mm)

l         = 0.2;
noX0      = 50;        % number of elements along X direction
noY0      = noX0;       % number of elements along Y direction (to have square cells)
ghostCell = 1;

[grid]    = buildGrid2D(l,l,noX0,noY0, ghostCell);

node      = grid.node;
element   = grid.element;
deltax    = grid.deltax;
deltay    = grid.deltay;
elemCount = grid.elemCount;
nodeCount = grid.nodeCount;
numx2     = grid.numx;
numy2     = grid.numy;


%% generate material points
ra            = 0.1/2;
ppc           = [3 3];
circle.center = [2*ra ra+grid.deltay];
circle.radius = ra;

[res]             = generateMPForCircle(circle,ppc,grid);
res.position(:,2) = res.position(:,2) + grid.deltay;

grid.lpx = deltax/ppc(1);
grid.lpy = deltay/ppc(2);

pCount        = size(res.position,1);
body1.volume  = res.volume;
body1.volume0 = res.volume;
body1.mass    = res.volume*density;
body1.coord   = res.position;
body1.coord0  = res.position;
body1.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress  = zeros(pCount,3);                % stress
body1.dStress = zeros(pCount,3);                % deviatoric stress
body1.strain  = zeros(pCount,3);                % strain
body1.velo    = zeros(pCount,2);                % velocity
body1.density = density*ones(pCount,1);         % updated density
body1.damage  = zeros(pCount,1);                % damage
body1.damThreshold = zeros(pCount,1);           % damage threshold
body1.eStrain = zeros(pCount,1);                % effective strain
body1.velo(:,2) = -v0;

% initialise minimum strain and damage
for ip=1:pCount
  Vol                    = body1.volume(ip);
  body1.damThreshold(ip) = (Vol*k)^(-1/m);
  
  xrand  = rand;
  %body1.damage(ip) = maxDamage0*(-log(rand))^(1/m) + minDamage0;
end

% put all bodies in one variable
bodies    = cell(1,1);
bodies{1} = body1;
bodyCount = length(bodies);

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

% find GIMP element connectivity

gimpElement = cell(elemCount,1);

for e=1:elemCount
    neighbors      = getNeighbors(e, grid.numx, grid.numy);
    neighborNodes  = element(neighbors,:);
    gimpElement{e} = unique(neighborNodes);
end

for ib=1:length(bodies)
  body      = bodies{ib};
  coord     = body.coord;
  elems     = ones(size(coord,1),1);
  
  for ip=1:size(coord,1)
    x = coord(ip,1); y = coord(ip,2);
    e = floor(x/deltax) + 1 + grid.numx*floor(y/deltay);
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


%% plot mesh, particles

coords=[bodies{:}.coord0];
hold on
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(bodies{1}.elements,:),'Q4','cy-',2.1);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','r-',2.1);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
%plot(res.position(:,1),res.position(:,2),'r.','markersize',10);
%plot(coords2(:,1),coords2(:,2),'r.','markersize',10);
axis off

%% node quantities

nmass      = zeros(nodeCount,1);  % nodal mass vector
nmomentum0 = zeros(nodeCount,2);  % nodal momentum vector
nmomentum  = zeros(nodeCount,2);  % nodal momentum vector
nmomentumS = zeros(nodeCount,2);  % nodal momentum vector, double mapping
niforce    = zeros(nodeCount,2);  % nodal internal force vector
neforce    = zeros(nodeCount,2);  % nodal external force vector



gBasis     = zeros(bodyCount,pCount,4);
gGradX     = zeros(bodyCount,pCount,4);
gGradY     = zeros(bodyCount,pCount,4);

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-14; % mass tolerance

c     = sqrt((bulkModulus)/density);
dtime = 0.1*grid.deltay/c;
time  = 10e-4;
t     = 0;

interval     = 2;
nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
  disp(['time step ',num2str(t)])
  
  %reset grid data (body contribution)
  nmass(:)     = 0;
  nmomentum0(:)= 0;
  nmomentumS(:)= 0;
  niforce(:)   = 0;
  neforce(:)   = 0;
  %% loop over bodies
  for ib=1:bodyCount
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints;
    for ie=1:length(elems)         % loop over computational cells or elements
      e     = elems(ie);
      esctr = gimpElement{e};
      enode = node(esctr,:);     % element node coords
      mpts  = mpoints{e};        % particles inside element e
      for p=1:length(mpts)       % loop over particles
        pid    = mpts(p);
        xp     = body.coord(pid,:);
        Mp     = body.mass(pid);
        vp     = body.velo(pid,:);
        Vp   = body.volume(pid);
        stress =  bodies{ib}.stress(pid,:);

        % loop over nodes of current element "ie"
        for i=1:length(esctr)
          id    = esctr(i);          
          x     = xp - node(id,:);
          [N,dNdx]=getGIMP2D(x,deltax,deltay,grid.lpx,grid.lpy);
          % store grid basis functions and gradients
          gBasis(ib,pid,i) = N;
          gGradX(ib,pid,i) = dNdx(:,1);
          gGradY(ib,pid,i) = dNdx(:,2);
        
          dNIdx = dNdx(1);
          dNIdy = dNdx(2);
          nmass(id)       = nmass(id)        + N*Mp;
          nmomentum0(id,:)= nmomentum0(id,:) + N*Mp*vp;
          niforce(id,1)   = niforce(id,1) - Vp*(stress(1)*dNIdx + stress(3)*dNIdy);
          niforce(id,2)   = niforce(id,2) - Vp*(stress(3)*dNIdx + stress(2)*dNIdy);
        end
      end
    end
  end
  
  % update momenta (including boundary conditions)
  nmomentum = nmomentum0 + dtime*niforce;
  nmomentum(grid.bNodes,:) = 0;
  
  %% update particle velocity and position
  %k = 0; u = 0;
  for ib=1:bodyCount
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints;
    % loop over computational cells or elements
    for ie=1:length(elems)
      e     = elems(ie);
      esctr = gimpElement{e};
      enode = node(esctr,:);     % element node coords
      mpts  = mpoints{e};       % particles inside element e
      % loop over particles
      for p=1:length(mpts)
        pid  = mpts(p);
        xp   = body.coord(pid,:);
        xp0  = body.coord0(pid,:);
        Mp   = body.mass(pid);
        vp   = body.velo(pid,:);
        Vp   = body.volume(pid);    
        for i=1:length(esctr)
          % retrieve the grid functions/grads
          %N    = gBasis(ib,pid,i);    
          id = esctr(i);
          x     = xp0 - node(id,:);
          [N,~]=getGIMP2D(x,deltax,deltay,grid.lpx,grid.lpy);
          if (nmass(id)~=0)
            massInv  = (1/nmass(id))*N;
            vp       = vp +  (nmomentum(id,:)-nmomentum0(id,:))*massInv;
            xp       = xp + dtime * nmomentum(id,:)*massInv;
          end
        end
        bodies{ib}.velo(pid,:) = vp;
        bodies{ib}.coord(pid,:)= xp;
        % mapped back bGrid momenta (used to compute L,,epsilon and stress)
        for i=1:length(esctr)
          id      = esctr(i);          
          x     = xp0 - node(id,:);
          [N,~]=getGIMP2D(x,deltax,deltay,grid.lpx,grid.lpy);
          nmomentumS(id,:)  = nmomentumS(id,:) +  Mp*vp*N;
        end
      end
    end
  end
  % do not forget boundary conditions on grid momenta
  nmomentumS(grid.bNodes,:) = 0;
  % update particle stresses
  for ib=1:bodyCount
    body      = bodies{ib};
    elems     = body.elements;
    mpoints   = body.mpoints;
    % loop over computational cells or elements
    for ie=1:length(elems)
      e     = elems(ie);
      esctr = gimpElement{e};
      enode = node(esctr,:);     % element node coords
      mpts  = mpoints{e};       % particles inside element e
      % loop over particles
      for p=1:length(mpts)
        pid  = mpts(p);
        xp0  = body.coord0(pid,:);        
        %vp   = body.velo(pid,:);        
        % retrieve the grid functions/grads
        %dNdx = [squeeze(gGradX(ib,pid,:)) squeeze(gGradY(ib,pid,:))];        
        Lp   = zeros(2,2);
        for i=1:length(esctr)
          id = esctr(i);          
          x     = xp0 - node(id,:);
          [~,dNdx]=getGIMP2D(x,deltax,deltay,grid.lpx,grid.lpy);
          if nmass(id) > 0
            vI  = nmomentumS(id,:)/nmass(id);  % nodal velocity
            Lp  = Lp + vI'*dNdx;
          end
        end
        % update stress last
        F       = ([1 0;0 1] + Lp*dtime)*reshape(bodies{ib}.deform(pid,:),2,2);
        bodies{ib}.deform(pid,:) = reshape(F,1,4);
        bodies{ib}.volume(pid  ) = det(F)*bodies{ib}.volume0(pid);
        if det(F) < 0, error('Negative Jacobian'); end
        dEps    = dtime * 0.5 * (Lp+Lp');
        traceEps= dEps(1,1) + dEps(2,2);
        dEpsDev = [dEps(1,1)-oneThird*traceEps dEps(2,2)-oneThird*traceEps dEps(1,2)];
        dShearStres = 2*bulkModulus*dEpsDev;
        newDensity  = bodies{ib}.density(pid)/(1.+traceEps);        
        bodies{ib}.density(pid) = newDensity;
        pressure   = c*c*(newDensity-density);
        dStress    = bodies{ib}.dStress(pid,:) + dShearStres;
        bodies{ib}.dStress(pid,:) = dStress;
        %dsigma      = -dPressure*I + dShearStres;
        % compute damage
        %sigma       = bodies{ib}.stress(pid,:)  + dsigma;        % undamaged stress (glob)
        sigma       = -pressure*I + dStress;
        bodies{ib}.stress(pid,:) = sigma;
        sigmaMatrix = [sigma(1) sigma(3);sigma(3) sigma(2)];     % to matrix format to use 'eig'
        [D,V]       = eig(sigmaMatrix);                          % comput principal stresses
        sigma1      = max(max(V));                               % maximum principal stress
        bodies{ib}.eStrain(pid) = (sigma1>0)*sigma1/(bulkModulus+4*shearModulus/3); % only tensile stress counts
        effStrain   = body.eStrain(pid);
        minStrain   = body.damThreshold(pid);                    % get damage threshold
        dam         = body.damage(pid);                          % get current damage
        if effStrain > minStrain
          %dam = dam + dtime*mAlpha*effStrain^mDiv3*dam^twoThird;
          p = [1 -dtime*mAlpha*effStrain^mDiv3 0 -dam];
          r = roots(p);
        end
        %dam=0;
        body.damage(pid) = dam;
        V(V>0) = (1-dam)*V(V>0);                  % (1-D)*sigma only for tensile stress
        sigma  = D'*V*D;                        % rotate the damaged stress back to global coord.
        bodies{ib}.stress(pid,:)  =  [sigma(1,1) sigma(2,2) sigma(1,2)];
      end
    end
  end
  
  % update the element particle list
  for ib=1:length(bodies)
    body      = bodies{ib};
    coord     = body.coord;
    elems     = ones(size(coord,1),1);
    
    for ip=1:size(coord,1)
      x = coord(ip,1); y = coord(ip,2);
      e = floor(x/deltax) + 1 + grid.numx*floor(y/deltay);
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
    bodies{ib}.coord0  = bodies{ib}.coord;
  end
  
  % store time,velocty for plotting
  %ta  = [ta;t];
  %xcm = [xcm;xx];
  
  % VTK output
  
  if (  mod(istep-1,interval) == 0 )
    xp = [bodies{:}.coord];
    s  = [bodies{:}.stress];
    stress = [s sum(s,2)/3];
    data.stress = stress;
    data.damage=[bodies{:}.damage];
    vtkFile = sprintf('%s%d',vtkFileName,istep-1);
    VTKParticles(xp,vtkFile,data);
  end
  
  
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

disp([num2str(toc),'   DONE '])
