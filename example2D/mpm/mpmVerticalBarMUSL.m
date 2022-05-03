% This file implements the Material Point Method of Sulsky 1994.
% Vertical bar with extreme deformation taken from
% A. Sadeghirad, R. M. Brannon, and J. Burghardt. A convected particle domain
% interpolation technique to extend applicability of the material point
% method for problems involving massive deformations. IJNME, 86(12):1435--1456, 2011.
%
% Modified Update Stress Last (MUSL) formulation.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% June 2014.

%%

addpath ../fem_util/
addpath ../fem-functions/
addpath ../post-processing/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%
E   = 1;           % Young's modulus
nu  = 0.3;         % Poisson ratio
rho = 1050e-12;    % density
K   = E/3/(1-2*nu);    % bulk modulus
mu    = E/2/(1+nu);% shear modulus
lambda = K - 2/3*mu;

g     = 1000e3; % gravity
bodyf = [0 -g];

I  = [1 0;0 1];

tic;

disp([num2str(toc),'   INITIALISATION '])

%%   particle distribution
%

l     = 1e3;
numx2 = 24;      % number of elements along X direction
numy2 = 24;      % number of elements along Y direction
[pmesh]= buildGrid2D(l,l,numx2,numy2, 0);

pmesh.node(:,1) = pmesh.node(:,1) +  7* l/2;
pmesh.node(:,2) = pmesh.node(:,2) + 25*l/2;

mpId = 12;% index of marked particle (either 3 or 12)

% store the particle mesh into a structure for convenience
elemType           = 'Q4';
particles.node     = pmesh.node;
particles.elem     = pmesh.element;
particles.elemType = elemType;

pCount  = pmesh.elemCount;                       % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
coords  = zeros(pCount,2);
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stressx
strain  = zeros(pCount,3);                % strain
velo    = zeros(pCount,2);

for e = 1:pmesh.elemCount
    coord = pmesh.node(pmesh.element(e,:),:);
    a     = 0.5*( coord(1,1)*coord(2,2)  - coord(2,1)*coord(1,2) ...
        + coord(2,1)*coord(3,2)  - coord(3,1)*coord(2,2) ...
        + coord(3,1)*coord(4,2)  - coord(4,1)*coord(3,2) ...
        + coord(4,1)*coord(1,2)  - coord(1,1)*coord(4,2) );
    volume(e)  = a;
    mass(e)    = a*rho;
    coords(e,:) = mean(coord); % center of each element=particle
end

volume0 = volume;
coords0 = coords;

lpx    = pmesh.deltax;
lpy    = pmesh.deltay;

dvec1   = zeros(pCount,2);                % domain vector 1, r1
dvec2   = zeros(pCount,2);                % domain vector 2, r2

for p=1:pCount
    dvec1(p,1) = 0.5*lpx;
    dvec1(p,2) = 0;
    dvec2(p,1) = 0;
    dvec2(p,2) = 0.5*lpy;
end

dvec10 = dvec1;
dvec20 = dvec2;

%% Computational grid

ghostCell=0;
lx     = (l/2)*16;
ly     = (l/2)*28;
numx2 = 16;      % number of elements along X direction
numy2 = 28;      % number of elements along Y direction
[mesh]= buildGrid2D(lx,ly,numx2,numy2, ghostCell);
element= mesh.element;
node   = mesh.node;


% find boundary nodes

fixNodes=find(abs(node(:,2)-27*l/2)<1e-10);

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

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum = zeros(mesh.nodeCount,2);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,2);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,2);  % nodal external force vector
nvelo     = zeros(mesh.nodeCount,2);  % nodal velocity vector (used for MUSL only)

%% plot mesh, particles

figure
hold on
plot_mesh(node,element,'Q4','k-',1.2);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(coords0(:,1),coords0(:,2),'k.','markersize',10);
plot(node(fixNodes,1),node(fixNodes,2),'r*','markersize',14);
axis off

ta = 0;
ka = 0;

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.1*(mesh.deltax/c);
time  = 0.25;
t     = 0;

istep = 1;

nsteps = floor(time/dtime);

pCoords    = cell(nsteps,1);

activeNodes = unique(element(unique(pElems),:));

while ( t < time )
    disp(['time step ',num2str(t)])
    % reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    neforce(:)   = 0;
    nvelo(:)     = 0;
    % loop over computational cells or elements
    for e=1:mesh.elemCount
        esctr = element(e,:);      % element connectivity
        enode = node(esctr,:);     % element node coords
        mpts  = mpoints{e};        % particles inside element e
        if (isempty(mpts)) continue;end
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            sigm = stress(pid,:);
            xp   = coords(pid,:);
            for i=1:length(esctr)
                id    = esctr(i);
                x       = xp - node(id,:);
                [N,dNdx]= getMPM2D(x,mesh.deltax,mesh.deltay);
                nmass(id)       = nmass(id)       + N*mass(pid);
                nmomentum(id,:) = nmomentum(id,:) + N*mass(pid)*velo(pid,:);
                niforce(id,1)   = niforce(id,1) - volume(pid)*(sigm(1)*dNdx(1) + sigm(3)*dNdx(2));
                niforce(id,2)   = niforce(id,2) - volume(pid)*(sigm(3)*dNdx(1) + sigm(2)*dNdx(2));
                neforce(id,:)   = neforce(id,:) + mass(p)*N*bodyf;
            end
        end
    end
    
    % update nodal momenta
    nforce    = niforce + neforce;
    nforce   (fixNodes,2)  = 0;
    nmomentum(fixNodes,2)  = 0;
    nmomentum = nmomentum + nforce*dtime;
    
    % update particle velocity
    for e=1:mesh.elemCount
        esctr = element(e,:);
        enode = node(esctr,:);
        mpts  = mpoints{e};
        if (isempty(mpts)) continue;end
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            xp   = coords(pid,:);
            for i=1:length(esctr)
                id = esctr(i);
                x       = xp - node(id,:);
                [N,dNdx]= getMPM2D(x,mesh.deltax,mesh.deltay);
                velo(pid,:)   = velo(pid,:) + dtime * N*nforce(id,:)  /nmass(id);
            end
            for i=1:length(esctr)
                id = esctr(i);
                x       = xp - node(id,:);
                [N,dNdx]= getMPM2D(x,mesh.deltax,mesh.deltay);
                nvelo(id,:)  = nvelo(id,:) + N*velo(pid,:)*mass(pid);
            end
        end
    end
    
    nvelo(activeNodes,1) = nvelo(activeNodes,1) ./ nmass(activeNodes);
    nvelo(activeNodes,2) = nvelo(activeNodes,2) ./ nmass(activeNodes);
    
    nvelo(fixNodes,2) = 0;
    
    % update particle position and stress
    for e=1:mesh.elemCount
        esctr = element(e,:);
        enode = node(esctr,:);
        mpts  = mpoints{e};
        if (isempty(mpts)) continue;end        
        for p=1:length(mpts)
            pid  = mpts(p);
            xp   = coords(pid,:);            
            Lp   = zeros(2,2);
            for i=1:length(esctr)
                id = esctr(i);
                vI = nvelo(id,:);
                x       = xp - node(id,:);
                [N,dNdx]= getMPM2D(x,mesh.deltax,mesh.deltay);
                coords(pid,:)  = coords(pid,:) + dtime * N*nmomentum(id,:)/nmass(id);
                Lp = Lp + vI'*dNdx;         % particle gradient velocity
            end
            
            F          = ([1 0;0 1] + Lp*dtime)*reshape(deform(pid,:),2,2);
            deform(pid,:)= reshape(F,1,4);
            volume(pid)  = det(F)*volume0(pid);
            J       = det(F);
            
            if (J<0) || ~isreal(J)
                error('negative or imagiary Jacobian J')
            end
            
            b       = F*F';
            sigma   = 1/J*( mu*(b-I) + lambda*log(J)*I );
            stress(pid,:)  = [sigma(1,1) sigma(2,2)  sigma(1,2) ];
        end
    end
    
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
    
    activeNodes = unique(element(unique(pElems),:));
    
    t     = t + dtime;
    istep = istep + 1;
    
    ta = [ta;t];
    ka = [ka;coords(mpId,2)-coords0(mpId,2)];
    
    pCoords{istep}=coords;
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])


cpdi=load('cpdiVerticalBar.mat');
mpm =load('mpm-VerticalBar.mat');

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(cpdi.ta(1:end),cpdi.ka(1:end),'r-','LineWidth',1.6);
plot(mpm.ta(1:end),mpm.ka(1:end),'black-','LineWidth',1.6);
xlabel('Time')
ylabel('Displacement')
legend('MPM/MUSL','CPDI','MPM/Cutoff')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
axis([0 0.25 0 -1800])

coords=pCoords{93};

figure
hold on
plot_mesh(node,element(e,:),'Q4','r-',1.2);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
%plot(coords0(:,1),coords0(:,2),'r.','markersize',10);
plot(node(fixNodes,1),node(fixNodes,2),'r*','markersize',14);
axis off

disp([num2str(toc),'   DONE '])
