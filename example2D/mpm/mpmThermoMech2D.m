% This file implements the  Material Point Method for thermo-elasticity.
% The grid: four-noded bilinear elements.
% Leapfrog time integration.
%.
%
% Vinh Phu Nguyen
% Monash University , Australia
% 1 October 2019.

%%
%parpool('open',8);

%%

addpath ../../grid/
addpath ../../basis/
addpath ../../particleGen/
addpath ../../constitutiveModels/
addpath ../../util/
addpath ../../geoMesh/
addpath ../../externals/
addpath ../../postProcessing/
addpath ../../mls/

%%
clc
clear
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
E      = 100;               % Young modulus
nu     = 0.2;               % Poisson ratio
rho    = 0.1;              % density
lambda = E*nu/(1+nu)/(1-2*nu);
mu     = E/2/(1+nu);
cv     = sqrt(E/rho);
k      = 10.;              % thermal conductivity
c      = 1.0;              % specific heat
alpha  = 0.25;             % thermal expansion coefficient
T0     = 0;                % initial temperature

option = 1;                % boundary particle option

identity  = [1 0;0 1];

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C           = elasticityMatrix(E,nu,stressState);
D           = inv(C);

vtkFileName  = 'mpm2DThermalMech';
interval     = 2;

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell = 0;
lx        = 5.8;
ly        = 5.8;
noX0      = 29;      % number of elements along X direction
noY0      = 29;      % number of elements along Y direction
[mesh]    = buildGrid2D(lx,ly,noX0,noY0, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;
omegaC    = mesh.deltax*mesh.deltay;

%% generate material points
ppc           = [2 2];
square.x      = [0 5]; % solid is a square 5 x 5
square.y      = [0 5];
[pmesh]       = buildGrid2D(5,5,25,25, ghostCell); % 25 cells <= (6/30) * 25 = 5 mm
[res]         = generateMPForRectangle(square,ppc,pmesh);

res.position(:,1) = res.position(:,1) + 2*mesh.deltax;
res.position(:,2) = res.position(:,2) + 2*mesh.deltay;

pCount  = size(res.position,1);
volume  = res.volume;
volume0 = res.volume;
mass    = res.volume*rho;
coords  = res.position;
coords0 = coords;
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stress
velo    = zeros(pCount,2);                % velocity
density = rho*ones(pCount,1);             % density
% thermal part
temp0   = zeros(pCount,1);                % temperature (old)
temp    = zeros(pCount,1);                % temperature
CC      = zeros(pCount,1);                % specifit heat
Q       = zeros(pCount,2);                % heat flux q

% initial temperature = 0, so do nothing
CC(:)   = c;

% initial velocities, initial stress=0
% for p=1:pCount
%   velo(p,1) = ...;
%   velo(p,2) = ...;
% end

% find nodes of Dirichlet boundary

% two options all particles in the boundary cell
% or just one column of particle

if ( option == 1 )
  idx = find(abs(coords(:,1)<3*mesh.deltax+0.001));
else
  idx = find(abs(coords(:,1)<2.5*mesh.deltax+0.001));
end

idx2 = intersect(find(abs(node(:,1)-(5+2*mesh.deltax)) < 1e-10),...
                 find(abs(node(:,2)-(5+2*mesh.deltay)) < 1e-10) );



idx1 = intersect(find(abs(node(:,1)-(5+2*mesh.deltax))<1e-10),...
                 find(abs(node(:,2)-(  2*mesh.deltay))<1e-10) );

%% plot mesh, particles
figure(1)
hold on
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(bodies{1}.elements,:),'Q4','cy-',2.1);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','r-',2.1);
plot(coords(:,1),coords(:,2),'r.','markersize',40);
plot(coords(idx,1),coords(idx,2),'blacks','markersize',20);
plot(node(idx1,1),node(idx1,2),'bs','markersize',20);
plot(node(idx2,1),node(idx2,2),'bs','markersize',20);
axis on

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

pElems  = ones(pCount,1);
mpoints = cell(mesh.elemCount,1);

for p=1:pCount
    x = coords(p,1);
    y = coords(p,2);
    e = floor(x/mesh.deltax) + 1 + mesh.numx*floor(y/mesh.deltay);
    pElems(p) = e;
end

for e=1:mesh.elemCount
    id  = find(pElems==e);
    mpoints{e}=id;
end


activeElems = unique(pElems);
activeNodes = unique(element(activeElems,:));

%% nodal quantities
nmass       = zeros(nodeCount,1);  % nodal mass vector
nmassT      = zeros(nodeCount,1);  % nodal mass vector for temperature
nmomentum   = zeros(nodeCount,2);  % nodal momentum vector (begin)
nmomentum0  = zeros(nodeCount,2);  % nodal momentum vector (final)
nmomentumT  = zeros(nodeCount,1);  % nodal momentum for temperature vector (begin)
nmomentumT0 = zeros(nodeCount,1);  % nodal momentum for temperature vector (final)
niforce     = zeros(nodeCount,2);  % nodal internal force vector
neforce     = zeros(nodeCount,2);  % nodal external force vector
niforceT    = zeros(nodeCount,1);  % nodal internal force vector due to temperature
% for MUSL
nmomentumMUSL   = zeros(nodeCount,2);  % nodal momentum vector (begin), mapped back
nmomentumTMUSL  = zeros(nodeCount,1);  % nodal momentum for temperature vector (begin)

nvelo   = zeros(nodeCount,2);  
nvelo0  = zeros(nodeCount,2);  
ntemp   = zeros(nodeCount,1);  
ntemp0  = zeros(nodeCount,1);  

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c     = sqrt(E/rho);
dtime = min(0.1*mesh.deltax/cv,0.1*mesh.deltax^2/2/k);
time  = 50*dtime;
t     = 0;

istep = 1;
nsteps = floor(time/dtime);
pDisp  = zeros(nsteps,1);
ta     = 0:dtime:time;

while ( t < time )
    disp(['time step ',num2str(t)])
    % reset grid data
    nmass(:)        = 0;
    nmassT(:)       = 0;
    niforce(:)      = 0;
    neforce(:)      = 0;
    niforceT(:)     = 0;
    nmomentum0(:)   = 0;
    nmomentumT0(:)  = 0;
    nmomentumMUSL(:)  = 0;
    nmomentumTMUSL(:) = 0;
    % loop over computational cells or elements
    for ie=1:length(activeElems)
        e     = activeElems(ie);
        esctr = element(e,:);
        enode = node(esctr,:);     % element node coords
        mpts  = mpoints{e};
        
        % external force and mass
        for p=1:length(mpts)
            pid  = mpts(p);
            xx   = coords(pid,:);       % position
            mm   = mass(pid);           % mass
            sigma = stress(pid,:);      % stress
            Vp    = volume(pid);        % volume
            vel   = velo(pid,:);        % velocity
            Tp    = temp(pid);          % temperature
            Cp    = CC(pid);            % specific heat
            Qp    = Q(pid,:);           % heat  flux
            for i=1:length(esctr)
                id              = esctr(i);
                x               = xx - node(id,:);
                [N,dNdx]        = getMPM2D(x,mesh.deltax,mesh.deltay);
                nmass(id)       = nmass(id)          + N*mm;
                nmassT(id)      = nmassT(id)         + N*mm*Cp;
                nmomentum0(id,:)= nmomentum0(id,:)   + N*mm*vel;
                nmomentumT0(id) = nmomentumT0(id)    + N*mm*Cp*Tp;
                % internal force due to stress
                niforce(id,1)   = niforce(id,1) - Vp*(sigma(1)*dNdx(1) + sigma(3)*dNdx(2));
                niforce(id,2)   = niforce(id,2) - Vp*(sigma(3)*dNdx(1) + sigma(2)*dNdx(2));
                % internal force due to temperature
                niforceT(id)    = niforceT(id) + Vp*dot(Qp,dNdx);
                %neforce(id,2)   = neforce(id,2) + mm*N*g;
            end
        end
        
    end
    
    % update nodal velocity
    nforce     = niforce + neforce;
    
    nmomentum  = nmomentum0  + dtime * nforce;
    nmomentumT = nmomentumT0 + dtime * niforceT;
    
    % calculate the grid velocity/temperature
    
    nvelo0(activeNodes,:) = nmomentum0(activeNodes,:) ./ nmass(activeNodes);
    nvelo(activeNodes,:)  = nmomentum(activeNodes,:)  ./ nmass(activeNodes);
    
    ntemp0(activeNodes)   = nmomentumT0(activeNodes) ./ nmassT(activeNodes);
    ntemp(activeNodes)    = nmomentumT(activeNodes)  ./ nmassT(activeNodes);
    
    % Dirichlet BCs for velopcity and temperature
    
    nvelo0(idx2,1) = 0; nvelo(idx2,1) = 0;
    nvelo0(idx1,2) = 0; nvelo(idx1,2) = 0;
    
    % enforcement of temperature Dirichlet on nodes not working
    %ntemp0(idx) = 1.; ntemp(idx) = 1.;
    
    
    % update particle velocity, temperature and map them back to the grid
    for ie=1:length(activeElems)
        e     = activeElems(ie);
        esctr = element(e,:);
        enode = node(esctr,:);
        mpts  = mpoints{e};
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            xp   = coords(pid,:);
            mm   = mass(pid);           % mass
            Cp   = CC(pid);             % specific heat
            for i=1:length(esctr)
                id           = esctr(i);
                x            = xp - node(id,:);
                [N,dNdx]     = getMPM2D(x,mesh.deltax,mesh.deltay);
                velo(pid,:)  = velo(pid,:) + N*(nvelo(id)-nvelo0(id) );
                temp(pid)    = temp(pid)   + N*(ntemp(id)-ntemp0(id) );
            end
            % mapping back to the grid
            for i=1:length(esctr)
                id           = esctr(i);
                x            = xp - node(id,:);
                [N,dNdx]     = getMPM2D(x,mesh.deltax,mesh.deltay);
                nmomentumMUSL(id,:)  = nmomentumMUSL(id,:)   + N*mm*velo(pid);
                nmomentumTMUSL(id)   = nmomentumTMUSL(id)  + N*mm*Cp*temp(pid,:);
            end
        end
    end
    
    % final updated grid velocity and temperature
    ntemp(activeNodes)    = nmomentumTMUSL(activeNodes)    ./ nmassT(activeNodes);
    nvelo(activeNodes,:)  = nmomentumMUSL(activeNodes,:)   ./ nmass(activeNodes);
    
    % Dirichlet BCs
    
    nvelo(idx2,1) = 0; nvelo(idx1,2) = 0;
    % enforcement of temperature Dirichlet on nodes not working
    %ntemp(idx) = 1.;
    temp(idx)   = 1.;
    
    % update particle position, stresses (G2P)
    for ie=1:length(activeElems)
        e     = activeElems(ie);
        esctr = element(e,:);
        enode = node(esctr,:);
        mpts  = mpoints{e};
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            xp   = coords(pid,:);
            Lp   = zeros(2,2); q =zeros(1,2);
            for i=1:length(esctr)
                id = esctr(i);
                %vI = [0 0];
                x             = xp - node(id,:);
                [N,dNdx]      = getMPM2D(x,mesh.deltax,mesh.deltay);
                coords(pid,:) = coords(pid,:) + dtime*N*nvelo(id,:);
                vI            = nvelo(id,:);           % nodal velocity
                Lp            = Lp + vI'*dNdx;         % particle gradient velocity
                
                % update particle flux q
                q = q - k * ntemp(id) * dNdx ;
            end
            
            % update particle flux q             
            Q(pid,:)      = q;
            
            F             = ([1 0;0 1] + Lp*dtime)*reshape(deform(pid,:),2,2);
            deform(pid,:) = reshape(F,1,4);
            volume(pid)   = det(F)*volume0(pid);
            %J             = det(F);
            %if (J<0), disp('error');end;
            %density(pid)  = rho/J;
            dEps          = dtime * 0.5 * (Lp+Lp');
            thermalStrain = alpha * ( temp(pid) - temp0(pid) ) * identity;
            elasticStrain = dEps - thermalStrain;
            dsigma        = C * [elasticStrain(1,1);elasticStrain(2,2);2*elasticStrain(1,2)] ;
            stress(pid,:) = stress(pid,:) + dsigma';

        end
    end
    
    temp0 = temp;
    
    
    % update the element particle list
    for p=1:pCount
        x = coords(p,1);
        y = coords(p,2);
        e = floor(x/mesh.deltax) + 1 + mesh.numx*floor(y/mesh.deltay);
        pElems(p) = e;
    end
    
    for e=1:mesh.elemCount
        id  = find(pElems==e);
        mpoints{e}=id;
    end
    
    activeElems = unique(pElems);
    activeNodes = unique(element(activeElems,:));
    
    t     = t + dtime;
    istep = istep + 1;
    
    if (  mod(istep-1,interval) == 0 )
        xp = coords;
       
        s1 = [stress zeros(size(stress,1),1)];
        data.stress  = [s1;];
        data.damage  = temp;
        vtkFile       = sprintf('../../results/mpm/thermal-mech2D/%s%d',vtkFileName,istep-1);
        VTKParticles(xp,vtkFile,data);
        pos{istep} = xp;
    end
    
end
%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

figure(1)
clf
plot_field(node,element,'Q4',ntemp);
hold on
%plot_mesh(node+scaleFact*[U(xs) U(ys)],element,elemType,'g.-',1);

disp([num2str(toc),'   DONE '])
