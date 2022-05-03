% This file implements the Generalized Interpolation Material Point Method
% 
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% USL with explicit time integration. Plane strain version.
% Rigid particles to emulate moving platen.
% Symmetric conditions on ghost nodes.
%
% Upsetting of a elastic-plastic billet.
% Example taken from Sulsky et al. 1994 paper.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% February 2014.
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

eye2x2 = [1 1 0; 
          1 1 0; 
          0 0 0];
I_dev  = eye(3) - 0.5*eye2x2;

%% Material properties
%

% billet
E   = 200e3;      % Young's modulus
nu  = 0.3;         % Poisson ratio
rho = 7800e-12;    % density
yield=700;          % yield stress
k1   = 300;           % Hardening modulus 
mu      = E/2/(1+nu);    % shear modulus
lambda  = E*nu/((1+nu)*(1-2*nu));
kappa   = lambda + mu;%bulk modulus

v0   = 100*1000;   %velocity of the platen

vtkFileName  = 'mpm2DUpSetting';
vtkFileName1 = '../results/gimp/upset/grid';
interval     = 5;


tic;

%% Computational grid

ghostCell= 1;
lx       = 24;
ly       = 18;
numx2    = 24;       % number of elements along X direction
numy2    = 18;      % number of elements along Y direction
[mesh]   = buildGrid2D(lx,ly,numx2,numy2, ghostCell);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;

%% generate material points
ppc    = 3;          % number of particles per cell one direction
l      = 10;
h      = 15;
numx2  = l*ppc;      % number of elements along X direction
numy2  = h*ppc;      % number of elements along Y direction
[pmesh]= buildGrid2D(l,h,numx2,numy2, 0);

pmesh.node(:,1) = pmesh.node(:,1) + ghostCell*mesh.deltax;
pmesh.node(:,2) = pmesh.node(:,2) + ghostCell*mesh.deltay;

pCount  = pmesh.elemCount;                       % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
coords  = zeros(pCount,2);

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

coords0 = coords;

body1.volume  = volume;
body1.volume0 = volume;
body1.mass    = mass;
body1.coord   = coords;
body1.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress  = zeros(pCount,3);                % stress
body1.strain  = zeros(pCount,3);                % strain
body1.pstrain = zeros(pCount,3);                % plastic strain
body1.alpha   = zeros(pCount,1);                % strain
body1.velo    = zeros(pCount,2);                % velocity

%2. rigid body that models the imposed velocity

l     = 24;
h     = mesh.deltay;
numx2 = l*ppc;      % number of elements along X direction
numy2 = 1*ppc;      % number of elements along Y direction
[rmesh]= buildGrid2D(l,h,numx2,numy2, 0);

rmesh.node(:,1) = rmesh.node(:,1) +      ghostCell*mesh.deltax;
rmesh.node(:,2) = rmesh.node(:,2) + 15 + ghostCell*mesh.deltay;

pCount  = rmesh.elemCount;                       % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
coords  = zeros(pCount,2);


for e = 1:rmesh.elemCount
    coord = rmesh.node(rmesh.element(e,:),:);
    a     = 0.5*( coord(1,1)*coord(2,2)  - coord(2,1)*coord(1,2) ...
        + coord(2,1)*coord(3,2)  - coord(3,1)*coord(2,2) ...
        + coord(3,1)*coord(4,2)  - coord(4,1)*coord(3,2) ...
        + coord(4,1)*coord(1,2)  - coord(1,1)*coord(4,2) );
    volume(e)  = a;
    mass(e)    = a*rho;
    coords(e,:) = mean(coord); % center of each element=particle
end


pCount = length(volume);

body2.volume = volume;
body2.volume0 = volume;
body2.mass   = mass;
body2.coord  = coords;
body2.deform = repmat([1 0 0 1],pCount,1);     % gradient deformation
body2.stress = zeros(pCount,3);                % stress
body2.strain = zeros(pCount,3);                % strain
body2.velo   = zeros(pCount,2);                % velocity
body2.velo(:,2) = -v0;

bodies    = cell(2,1);
bodies{1} = body1;
bodies{2} = body2;
bodyCount = length(bodies);

% particle size (constant in uGIMP)
lpx    = pmesh.deltax;
lpy    = pmesh.deltay;

%% find GIMP element connectivity

gimpElement = cell(elemCount,1);

for e=1:elemCount
  neighbors      = getNeighbors(e, mesh.numx, mesh.numy);
  neighborNodes  = element(neighbors,:);
  gimpElement{e} = unique(neighborNodes);
end

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

for ib=1:length(bodies)
    body      = bodies{ib};
    coord     = body.coord;
    elems     = ones(size(coord,1),1);
    
    for ip=1:size(coord,1)
        x = coord(ip,1); y = coord(ip,2);
        e = floor(x/mesh.deltax) + 1 + mesh.numx*floor(y/mesh.deltay);
        elems(ip) = e;
    end
    
    bodies{ib}.elements = unique(elems);
    bodies{ib}.nodes    = unique(vertcat(gimpElement{bodies{ib}.elements}));
    
    if ib==2
      bodies{ib}.nodes    = unique(element(bodies{ib}.elements,:));
    end
    
    mpoints = cell(elemCount,1);
    for ie=1:elemCount
        id  = find(elems==ie);
        mpoints{ie}=id;
    end
    
    bodies{ib}.mpoints  = mpoints;  
end


    
%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,2);  % nodal momentum vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector
nvelo     = zeros(nodeCount,2);  % nodal velocities

%% find boundary conditions

eps=1e-12;
bottom0 = find(abs(node(:,2))<eps);
bottom1 = find(abs(node(:,2)-mesh.deltay)<eps);
bottom2 = find(abs(node(:,2)-2*mesh.deltay)<eps);
left0   = find(abs(node(:,1))<eps);
left1   = find(abs(node(:,1)-mesh.deltax)<eps);
left2   = find(abs(node(:,1)-2*mesh.deltax)<eps);

fixedY    = [bottom1; bottom2];
fixedX    = [left1; left2];

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
plot(node(bottom0,1),node(bottom0,2),'cy.','markersize',15);
plot(node(bottom1,1),node(bottom1,2),'cy*','markersize',15);
plot(node(bottom2,1),node(bottom2,2),'cys','markersize',15);
plot(node(left0,1),node(left0,2),'r.','markersize',15);
plot(node(left1,1),node(left1,2),'r*','markersize',15);
plot(node(left2,1),node(left2,2),'rs','markersize',15);
%plot(node(bodies{1}.nodes,1),node(bodies{1}.nodes,2),'b.','markersize',15);
%axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 5e-8;
time  = 10e-5;
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
    disp(['time step ',num2str(t)])
    %% reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    %nvelo(:)     = 0;        
    for ib=1:1                         %% loop over bodies 
        body      = bodies{ib};
        elems     = body.elements;
        mpoints   = body.mpoints;        
        for ie=1:length(elems)         % loop over computational cells or elements
            e     = elems(ie);
            esctr = gimpElement{e};    % element connectivity            
            mpts  = mpoints{e};        % particles inside element e                             
            for p=1:length(mpts)       % loop over particles 
                pid    = mpts(p);
                xp     = body.coord(pid,:);
                stress = body.stress(pid,:);
                Mp     = body.mass(pid);
                vp     = body.velo(pid,:);
                Vp     = body.volume(pid);   
                for i=1:length(esctr)  % loop over nodes of current element "ie" 
                    id    = esctr(i);
                    x     = xp - node(id,:);
                    [N,dNdx]=getGIMP2D(x,mesh.deltax,mesh.deltay,lpx,lpy);                    
                    dNIdx = dNdx(1);
                    dNIdy = dNdx(2);
                    nmass(id)       = nmass(id)       + N*Mp;
                    nmomentum(id,:) = nmomentum(id,:) + N*Mp*vp;
                    niforce(id,1)   = niforce(id,1) - Vp*(stress(1)*dNIdx + stress(3)*dNIdy);
                    niforce(id,2)   = niforce(id,2) - Vp*(stress(3)*dNIdx + stress(2)*dNIdy);                                        
                end
            end
        end
        
        % update nodal momenta
        
        activeNodes=bodies{ib}.nodes;
        nmomentum(activeNodes,:) = nmomentum(activeNodes,:) + niforce(activeNodes,:)*dtime;
    end
        
    % boundary conditions at symmetric lines                  
    nmomentum(bottom1,2) = 0; % fixed boundary conditions
    nmomentum(left1,1)   = 0; % fixed boundary conditions       
    niforce(left1,1)     = 0; % fixed boundary conditions       
    niforce(bottom1,2)   = 0; % fixed boundary conditions   
    
    % symmetric conditions at ghost nodes
    nmass(bottom0) = nmass(bottom2);
    nmass(left0)   = nmass(left2);
    
    nmomentum(bottom0,2) = -nmomentum(bottom2,2);
    nmomentum(left0,1)   = -nmomentum(left2,1);
    
    niforce(bottom0,2) = -niforce(bottom2,2);
    niforce(left0,1)   = -niforce(left2,1);
    
    % boundary conditions from rigid particles
    %nmomentum(bodies{2}.nodes,1) = 0;
    nmomentum(bodies{2}.nodes,2) = nmass(bodies{2}.nodes)*(-v0);
    niforce(bodies{2}.nodes,2)   = 0; 
    
    for ib=1:1
        body      = bodies{ib};
        elems     = body.elements;
        mpoints   = body.mpoints;        
        for ie=1:length(elems)         % loop over computational cells or elements
            e     = elems(ie);
            esctr = gimpElement{e};    % element connectivity            
            mpts  = mpoints{e};        % particles inside element e            
            for p=1:length(mpts)       % loop over particles
                pid  = mpts(p);
                xp   = body.coord(pid,:);
                xp0  = body.coord(pid,:);                                
                vp   = body.velo(pid,:);                
                Lp   = zeros(2,2);
                for i=1:length(esctr)
                    id  = esctr(i);
                    x   = xp0 - node(id,:);
                    [N,dNdx]=getGIMP2D(x,mesh.deltax,mesh.deltay,lpx,lpy);
                    if (nmass(id)~=0)
                      xp  = xp + dtime * N*nmomentum(id,:)/nmass(id);
                      vp  = vp + dtime * N*niforce(id,:)/nmass(id);
                      Lp  = Lp + (nmomentum(id,:)/nmass(id))'*dNdx;% particle gradient velocity
                    end
                end
                                
                bodies{ib}.coord(pid,:)= xp;
                bodies{ib}.velo(pid,:) = vp;
                
                F       = ([1 0;0 1] + Lp*dtime)*reshape(bodies{ib}.deform(pid,:),2,2);
                bodies{ib}.deform(pid,:) = reshape(F,1,4);
                bodies{ib}.volume(pid  ) = det(F)*bodies{ib}.volume0(pid);
                dEps    = dtime * 0.5 * (Lp+Lp');
                bodies{ib}.strain(pid,:)  = bodies{ib}.strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
                
                if ( ib == 2 ) % elastic
%                     dsigma  = C2 * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
%                     bodies{ib}.stress(pid,:)  = bodies{ib}.stress(pid,:) + dsigma';
                else           % plastic
                     
                    epsilon   =  bodies{ib}.strain (pid,:)';   % column vector
                    epsilonp0 =  bodies{ib}.pstrain(pid,:)';
                    alpha0    =  bodies{ib}.alpha  (pid);
                    
                    [sigma,alpha,epsilonp] = updateVonMisesMaterial ( epsilon,epsilonp0,alpha0,mu,kappa,k1,yield );
                    
                    bodies{ib}.stress (pid,:) = sigma;
                    bodies{ib}.pstrain(pid,:) = epsilonp;
                    bodies{ib}.alpha  (pid,:) = alpha;
                end                                                
            end
        end
    end
    
    bodies{2}.coord = bodies{2}.coord + dtime* bodies{2}.velo;
    
    % update the element particle list
    for ib=1:length(bodies)
        body      = bodies{ib};
        coord     = body.coord;
        elems     = ones(size(coord,1),1);
        
        for ip=1:length(elems)
            x = coord(ip,1); y = coord(ip,2);
            e = floor(x/mesh.deltax) + 1 + mesh.numx*floor(y/mesh.deltay);
            elems(ip) = e;
        end
        
        bodies{ib}.elements = unique(elems);
        bodies{ib}.nodes    = unique(vertcat(gimpElement{bodies{ib}.elements}));
%         if ib==2
%           bodies{ib}.nodes    = unique(element(bodies{ib}.elements,:));
%         end
        mpoints = cell(elemCount,1);
        for ie=1:elemCount
            id  = find(elems==ie);
            mpoints{ie}=id;
        end
        
        bodies{ib}.mpoints  = mpoints;
    end
        
    % VTK output    
    if (  mod(istep-1,interval) == 0 )
        xp = [bodies{1}.coord;bodies{2}.coord];
        s2 = [bodies{2}.stress zeros(size(bodies{2}.stress,1),1)];
        s1 = [bodies{1}.stress zeros(size(bodies{1}.stress,1),1)];
        for i=1:size(bodies{1}.stress,1)
            devSig = I_dev*s1(i,1:3)';
            s1(i,4)= sqrt(  devSig(1)^2 + devSig(2)^2 + 2*devSig(3)^2  );
        end
        data.stress  = [s1;s2];
        data.pstrain = [bodies{1}.alpha;zeros(size(bodies{2}.stress,1),1)];
        data.velo    = [bodies{1}.velo; bodies{1}.velo];
        vtkFile = sprintf('../results/gimp/upset/%s%d',vtkFileName,istep-1);
        VTKParticles(xp,vtkFile,data);
        pos{istep} = xp;
    end
        
    % advance to the next time step    
    t = t + dtime;
    istep = istep + 1;        
end

Ux= zeros(size(mesh.node,1),1);
Uy= zeros(size(mesh.node,1),1);
sigmaXX = zeros(size(mesh.node,1),1);
sigmaYY = zeros(size(mesh.node,1),1);
sigmaXY = zeros(size(mesh.node,1),1);

VTKPostProcess(mesh.node,mesh.element,2,'Quad4',vtkFileName1,...
    [sigmaXX sigmaYY sigmaXY],[Ux Uy]);
  
disp([num2str(toc),'   DONE '])
