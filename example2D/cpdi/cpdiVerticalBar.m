% This file implements the Material Point Method with CPDI interpolation
% described in the article
%
% A. Sadeghirad, R. M. Brannon, and J. Burghardt. A convected particle domain
% interpolation technique to extend applicability of the material point
% method for problems involving massive deformations. IJNME, 86(12):1435--1456, 2011.

% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Particles taken from a Q4 unstructured mesh (Gmsh).
%
% Vertical bar with extreme deformation.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% May 2014, Saigon, Vietnam.

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

g     = 1500e3; % gravity
bodyf = [0 -g];

I  = [1 0;0 1];

interval     = 100;
vtkFileName  = 'bar';
vtkFileName1 = '../results/cpdi/verticalBar/barGrid';

tic;

disp([num2str(toc),'   INITIALISATION '])

%%   particle distribution
%

l     = 1e3;
numx2 = 6;      % number of elements along X direction
numy2 = 6;      % number of elements along Y direction
[pmesh]= buildGrid2D(l,l,numx2,numy2, 0);

pmesh.node(:,1) = pmesh.node(:,1) +   l/2;
pmesh.node(:,2) = pmesh.node(:,2) + 18*l/2;

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
lx     = (l/2)*4;
ly     = (l/2)*21;
numx2 = 4;      % number of elements along X direction
numy2 = 21;      % number of elements along Y direction
[mesh]= buildGrid2D(lx,ly,numx2,numy2, ghostCell);
element= mesh.element;
node   = mesh.node;


% find boundary nodes

fixNodes=find(abs(node(:,2)-20*l/2)<1e-10);

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum = zeros(mesh.nodeCount,2);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,2);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,2);  % nodal external force vector
nforce    = zeros(mesh.nodeCount,2);  % nodal force vector

%% plot mesh, particles

hold on
plot_mesh(particles.node,particles.elem,elemType,'b-',0.2);
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
plot(node(fixNodes,1),node(fixNodes,2),'r*','markersize',14);
axis off

ta = 0;           % time
ka = 0;           % kinetic energy


%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-12; % mass tolerance

c     = sqrt(E/rho);
% for g=3500 use CFL=0.01
% for g=1000 use CFL=0.1
dtime = 0.1*(mesh.deltax/c);
time  = 0.25;
t     = 0;


nsteps = floor(time/dtime);

pmeshOut   = cell(nsteps,1);
pCoords    = cell(nsteps,1);

istep = 1;

while ( t < time )
    disp(['time step ',num2str(t)])
    % reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    neforce(:)   = 0;
    % loop over particles
    for p=1:pCount
        sig    = stress(p,:);
        xp     = coords(p,:);
        r1p    = dvec1(p,:);
        r2p    = dvec2(p,:);
        % particle mass and momentum to node
        data  = cpdi2D(xp,r1p,r2p,mesh);
        esctr = data.node;
        for i=1:length(esctr)
            id              = esctr(i);
            nmass(id)       = nmass(id)       + data.phi(i)*mass(p);
            nmomentum(id,:) = nmomentum(id,:) + data.phi(i)*mass(p)*velo(p,:);
            niforce(id,1)   = niforce(id,1) - volume(p)*(sig(1)*data.dphi(i,1) + sig(3)*data.dphi(i,2));
            niforce(id,2)   = niforce(id,2) - volume(p)*(sig(3)*data.dphi(i,1) + sig(2)*data.dphi(i,2));
            %if (t==0)
            neforce(id,:)   = neforce(id,:) + mass(p)*data.phi(i)*bodyf;
            %end
        end
    end
    
    % update nodal momenta
    nforce    = niforce + neforce;
    nforce   (fixNodes,2)  = 0;
    nmomentum(fixNodes,2)  = 0;
    
    nmomentum = nmomentum + nforce*dtime;
    
    % update particle velocity and position and stresses
    
    % loop over particles
    for p=1:pCount
        Lp    = zeros(2,2);
        xp    = coords(p,:);
        r1p   = dvec1(p,:);
        r2p   = dvec2(p,:);
        data  = cpdi2D(xp,r1p,r2p,mesh);
        esctr = data.node;
        for i=1:length(esctr)
            id = esctr(i);
            vI = [0 0];
            if nmass(id) > tol
                velo(p,:)  = velo(p,:)  + dtime * data.phi(i)*nforce(id,:) /nmass(id);
                coords(p,:)= coords(p,:)+ dtime * data.phi(i)*nmomentum(id,:)/nmass(id);
                vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
            end
            Lp = Lp + vI'*data.dphi(i,:);         % particle gradient velocity
        end
        
        F          = ([1 0;0 1] + Lp*dtime)*reshape(deform(p,:),2,2);
        deform(p,:)= reshape(F,1,4);
        dvec1(p,:) = dvec10(p,:)*F';
        dvec2(p,:) = dvec20(p,:)*F';
        volume(p)  = det(F)*volume0(p);
        J       = det(F);
        b       = F*F';
        sigma   = 1/J*( mu*(b-I) + lambda*log(J)*I );
        stress(p,:)  = [sigma(1,1) sigma(2,2)  sigma(1,2) ];
    end
    
    % store time,velocity for plotting
    
    
    % VTK output
    
%     if (  mod(istep,interval) == 0 )
%         vtkFile = sprintf('../results/cpdi/verticalBar/%s%d',vtkFileName,istep);
%         data.stress  = [stress zeros(pCount,1)];
%         data.pstrain = [];
%         data.velo    = velo;
%         VTKParticlesCPDI(particles,vtkFile,data);
%     end
    
    % write particle domains
    xx=[];
    for p=1:pCount
        xp    = coords(p,:);
        r1p   = dvec1(p,:);
        r2p   = dvec2(p,:);
        
        x1 = xp - r1p - r2p;
        x2 = xp + r1p - r2p;
        x3 = xp + r1p + r2p;
        x4 = xp - r1p + r2p;
        
        xx = [xx;x1;x2;x3;x4];                        
    end
    
    pmeshOut{istep} = xx;
    pCoords{istep}  = coords;
    
    % advance to the next time step
    
    t     = t + dtime;
    istep = istep + 1;
    
    ta = [ta;t];
    ka = [ka;coords(3,2)-coords0(3,2)];
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

% pvdFile = fopen(strcat('../results/cpdi2/',vtkFileName,'.pvd'), 'wt');
% 
% fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
% fprintf(pvdFile,'<Collection>\n');
% 
% for i = 1:nsteps
%     if (  mod(i,interval) == 0 )
%         vtuFile = sprintf('%s%d%s',vtkFileName,i,'.vtu');
%         fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
%     end
% end
% 
% fprintf(pvdFile,'</Collection>\n');
% fprintf(pvdFile,'</VTKFile>\n');
% 
% fclose(pvdFile);
% 
% %%
% Ux= zeros(size(mesh.node,1),1);
% Uy= zeros(size(mesh.node,1),1);
% sigmaXX = zeros(size(mesh.node,1),1);
% sigmaYY = zeros(size(mesh.node,1),1);
% sigmaXY = zeros(size(mesh.node,1),1);
% 
% VTKPostProcess(mesh.node,mesh.element,2,'Quad4',vtkFileName1,...
%     [sigmaXX sigmaYY sigmaXY],[Ux Uy]);

%%

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end)/1000,'b-','LineWidth',1.6);
xlabel('time')
ylabel('displacement [m]')
%legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
axis([0 0.25 0 -9])

save('cpdiVerticalBar.mat','ta','ka');

%%

step=60;

corners = pmeshOut{step};
coords  = pCoords {step};

figure
hold on
plot_mesh(node,element,'Q4','k-',1.2);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
for p=1:pCount
    plot_mesh(corners(4*(p-1)+1:4*p,:),[1 2 3 4],'Q4','r-',1.);
end
axis off

disp([num2str(toc),'   DONE '])



