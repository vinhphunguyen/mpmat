% This file implements the Material Point Method with CPDI interpolation
% described in the article
%
% A. Sadeghirad, R. M. Brannon, and J. Burghardt. A convected particle domain
% interpolation technique to extend applicability of the material point
% method for problems involving massive deformations. IJNME, 86(12):1435--1456, 2011.

% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Particles taken from a Q4 structured mesh.
%
% Vibration of a compilant cantilever beam
%
% Vinh Phu Nguyen
% Monash University
% 23 December 2016.

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
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties

E0   = 1e6;            % Young's modulus
nu0  = 0.3;            % Poisson ratio
rho0 = 1050;            % density
kappa = 3-4*nu0;       % Kolosov constant
mu    = E0/2/(1+nu0);  % shear modulus
K   = E0/3/(1-2*nu0);    % bulk modulus
lambda = K - 2/3*mu;

g     = 10;          % gravity m/s^2

I = [1 0;0 1];

bodyf = [0 -g];


interval  = 5;
vtkFileName  = 'cpdiQ4CantileverBeam';
vtkFileName1  = '../../results/cpdi/cpdiCantileverBeamGrid';


tic;

disp([num2str(toc),'   INITIALISATION '])



%% Computational grid

ghostCell=0;
lx     = 8;
ly     = 8;
numx2 = 16;      % number of elements along X direction
numy2 = 16;      % number of elements along Y direction
[mesh]= buildGrid2D(lx,ly,numx2,numy2, ghostCell);
element= mesh.element;
node   = mesh.node;


% find boundary nodes

fixNodes=find(abs(node(:,1)-mesh.deltax)<1e-10);

%%   particle distribution
%

l     = 4;
h     = 1;
numx2 = 24;      % number of elements along X direction
numy2 = 6;       % number of elements along Y direction
[pmesh]= buildGrid2D(l,h,numx2,numy2, 0);

% move the solid to the correct position
pmesh.node(:,1) = pmesh.node(:,1) +   mesh.deltax;
pmesh.node(:,2) = pmesh.node(:,2) + 6;

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
color   = ones(pCount,1);

for e = 1:pmesh.elemCount
    coord = pmesh.node(pmesh.element(e,:),:);
    a     = 0.5*( coord(1,1)*coord(2,2)  - coord(2,1)*coord(1,2) ...
                + coord(2,1)*coord(3,2)  - coord(3,1)*coord(2,2) ...
                + coord(3,1)*coord(4,2)  - coord(4,1)*coord(3,2) ...
                + coord(4,1)*coord(1,2)  - coord(1,1)*coord(4,2) );
    volume(e)  = a;
    mass(e)    = a*rho0;
    coords(e,:) = mean(coord); % center of each element=particle  
end

volume0 = volume;

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum = zeros(mesh.nodeCount,2);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,2);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,2);  % nodal external force vector
nforce    = zeros(mesh.nodeCount,2);  % nodal force vector

%% plot mesh, particles

hold on
plot_mesh(particles.node,particles.elem,elemType,'black-',1.2);
plot_mesh(node,element,'Q4','r-',1.2);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
plot(node(fixNodes,1),node(fixNodes,2),'r*','markersize',14);
axis off

ta = 0;           % time
ka = 0;

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c     = sqrt(E0/rho0);
dtime = 0.2*(mesh.deltax/c);
time  = 3.;
t     = 0;
nsteps = floor(time/dtime);

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
        % particle mass and momentum to node
        data  = cpdi22D(p,particles,mesh);
        esctr = data.node;
        for i=1:length(esctr)
            id              = esctr(i);
            nmass(id)       = nmass(id)       + data.phi(i)*mass(p);
            nmomentum(id,:) = nmomentum(id,:) + data.phi(i)*mass(p)*velo(p,:);
            niforce(id,1)   = niforce(id,1) - volume(p)*(sig(1)*data.dphi(i,1) + sig(3)*data.dphi(i,2));
            niforce(id,2)   = niforce(id,2) - volume(p)*(sig(3)*data.dphi(i,1) + sig(2)*data.dphi(i,2));            
            neforce(id,:)   = neforce(id,:) + mass(p)*data.phi(i)*bodyf;           
        end
    end
        
    % update nodal momenta   
    nforce    = niforce + neforce;
    nforce   (fixNodes,1)  = 0;
    nforce   (fixNodes,2)  = 0;
    nmomentum(fixNodes,1)  = 0;
    nmomentum(fixNodes,2)  = 0;
    
    nmomentum = nmomentum + nforce*dtime;
        
    % update particle velocity and position and stresses
    
    % loop over particles
    for p=1:pCount
        Lp    = zeros(2,2);        
        data  = cpdi22D(p,particles,mesh);
        esctr = data.node;
        for i=1:length(esctr)
            id = esctr(i);
            vI = [0 0];
            if nmass(id) > tol
                velo(p,:)  = velo(p,:) + dtime * data.phi(i)*nforce(id,:)/nmass(id);                
                vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
            end
            Lp = Lp + vI'*data.dphi(i,:);         % particle gradient velocity
        end
        
        F          = ([1 0;0 1] + Lp*dtime)*reshape(deform(p,:),2,2);
        deform(p,:)= reshape(F,1,4);
        volume(p)  = det(F)*volume0(p);
        J       = det(F);
        b       = F*F';
        sigma   = 1/J*( mu*(b-I) + lambda*log(J)*I );
        stress(p,:)  = [sigma(1,1) sigma(2,2)  sigma(1,2) ];
    end
    
    for c=1:size(particles.node,1)        
        xc    = particles.node(c,:);
        ec    = point2ElemIndex(xc,mesh);
        esctr = element(ec,:);
        for i=1:length(esctr)
            id      = esctr(i);
            x       = xc - node(id,:);
            [N,dNdx]= getMPM2D(x,mesh.deltax,mesh.deltay);
            if nmass(id) > tol                
                xc = xc + dtime*N*nmomentum(id,:)/nmass(id);                
            end            
        end
        particles.node(c,:) = xc;
    end
    
  
    % VTK output
    
    if (  mod(istep,interval) == 0 )
        vtkFile = sprintf('../../results/cpdi/%s%d',vtkFileName,istep);
        data.stress  = [stress];
        data.pstrain = [];
        data.color   = color;
        data.velo    = velo;
        VTKParticlesCPDI(particles,vtkFile,data);
    end
    
    
    % advance to the next time step
    
    t     = t + dtime;
    istep = istep + 1;
    
    % store time,velocty for plotting
    coord   = particles.node(particles.elem(numx2,:),:);
    coords3 = mean(coord); % center of each element=particle
    u  = coords3(:,2)-coords(numx2,2);
    ta = [ta;t];
    ka = [ka;u];        
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

pvdFile = fopen(strcat('../../results/cpdi/',vtkFileName,'.pvd'), 'wt');

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');

for i = 1:nsteps
    if (  mod(i,interval) == 0 )
        vtuFile = sprintf('%s%d%s',vtkFileName,i,'.vtu');
        fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
    end
end

fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);

%%
Ux= zeros(size(mesh.node,1),1);
Uy= zeros(size(mesh.node,1),1);
sigmaXX = zeros(size(mesh.node,1),1);
sigmaYY = zeros(size(mesh.node,1),1);
sigmaXY = zeros(size(mesh.node,1),1);

VTKPostProcess(mesh.node,mesh.element,2,'Quad4',vtkFileName1,...
    [sigmaXX sigmaYY sigmaXY],[Ux Uy]);

%%



ss=load('cpdi2VerticalBar.mat');

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b*','LineWidth',1.6);
plot(ss.ta(1:end),ss.ka(1:end),'r-','LineWidth',1.6);
xlabel('Time')
ylabel('Displacement')
legend('CPDI2','CPDI')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
axis([0 0.25 0 -1800])

save('cpdi2VerticalBar.mat', 'ta', 'ka');
%%
for e = 1:pmesh.elemCount
    coord = particles.node(particles.elem(e,:),:);
    coords(e,:) = mean(coord); % center of each element=particle
end

figure
hold on
plot_mesh(particles.node,particles.elem,elemType,'r-',0.2);
plot_mesh(node,element,'Q4','k-',1.)
plot(coords(:,1),coords(:,2),'k.','markersize',10);
axis off

disp([num2str(toc),'   DONE '])



