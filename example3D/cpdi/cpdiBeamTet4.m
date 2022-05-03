% This file implements the Material Point Method with CPDI2 interpolation
% described in the article
%
% A. Sadeghirad, R. M. Brannon, and J.E. Guilkey. Second-order convected
% particle domain in- terpolation (CPDI2) with enrichment for weak discontinuities
% at material interfaces. IJNME, 95(11):928-952, 2013.
%
% Three dimensional problems.
% The grid is a structured mesh consisting of 8-noded trilinear elements (H8).
% CPDI-Tet4: linear tetrahedron
%
% Compliant cantilever beam with extreme deformation.
%
% Vinh Phu Nguyen
% Monash University
% 13 April, 2016 (my son's visa granted this day)

%%
addpath ../../nurbs/nurbs-geopdes/inst
addpath ../../nurbs/nurbs-util/
addpath ../../util/
addpath ../../fem/
addpath ../../postProcessing/
addpath ../../constitutiveModels/
addpath ../../grid/
addpath ../../basis/
addpath ../../geoMesh/


%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

%Unit: m/s
%% Material properties
%

E     = 1e6;           % Young's modulus
nu    = 0.3;           % Poisson ratio
rho   = 1050;          % density
K     = E/3/(1-2*nu);  % bulk modulus
mu    = E/2/(1+nu);  % shear modulus
lambda = K - 2/3*mu;

g     = 10;          % gravity
bodyf = [0 -g 0];

I = [1 0 0;0 1 0; 0 0 1]; % identity matrix

% be careful with vtkFileName1 and change it according to your computer!!!
interval     = 10;
vtkFileName  = 'cantilever';
vtkFileName1 = '../../results/cpdi/cantileverGrid';

tic;

disp([num2str(toc),'   INITIALISATION '])

%%   particle distribution from a mesh
%

%meshFile = 'beam1000.msh'; % 1000 elements, works
meshFile = 'beam600.msh'; %491 tetrahedra, also works
%meshFile = 'beam178.msh'; %178 tetrahedra, also works, less accurate though
%meshFile = 'beam82.msh'; %82 tetrahedra, also works, less accurate though
gmesh    = load_gmsh (meshFile);

elemType = 'H4';
numnode  = gmesh.nbNod;
numelem  = gmesh.nbTets;
node1    = gmesh.POS(:,1:3);
element1 = gmesh.TETS(1:numelem,1:4);

% marked corner (recording displacement)

markedNode1 = find(abs(node1(:,1)-4)<1e-10);
markedNode2 = find(abs(node1(:,2)-0)<1e-10);
markedNode3 = find(abs(node1(:,3)-0)<1e-10);

markedNode  = intersect(intersect(markedNode1,markedNode2),markedNode3);

pCount  = numelem;                       % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
coords  = zeros(pCount,3);
deform  = repmat([1 0 0 0 1 0 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,6);                % stress
strain  = zeros(pCount,6);                % strain
velo    = zeros(pCount,3);
color   = ones (pCount,1);

[W,Q]   = quadrature( 1, 'TRIANGULAR', 3 ); 

for e = 1:numelem
    corners = node1(element1(e,:),:);
    t = fliplr(corners');
    corners=t';
    x1   = corners(1,1); y1   = corners(1,2); z1   = corners(1,3);
    x2   = corners(2,1); y2   = corners(2,2); z2   = corners(2,3);
    x3   = corners(3,1); y3   = corners(3,2); z3   = corners(3,3);
    x4   = corners(4,1); y4   = corners(4,2); z4   = corners(4,3);
    
    x21 = x2 - x1; x32 = x3 - x2; x43 = x4 - x3; x42 = x4 - x2;
    y23 = y2 - y3; y34 = y3 - y4; y12 = y1 - y2; y42 = y4 - y2;
    z34 = z3 - z4; z23 = z2 - z3; z12 = z1 - z2; z42 = z4 - z2;
    
    V   = (1/6)*(x21*(y23*z34-y34*z23) + x32*(y34*z12-y12*z34) + x43*(y12*z23-y23*z12));
%     V
%     
%     [N,dNdxi]= lagrange_basis('H4',Q(1,:));   % element shape functions    
%     J0       = corners'*dNdxi;                  % element Jacobian matrix
%     detJ     = det(J0)*W(1)
    
    volume(e)  = V;
    mass(e)    = V*rho;
    %coords(e,:) = mean(coord); % center of each element=particle  
end

volume0 = volume;

%% Computational grid

ghostCell=0;
lx     = 8;
ly     = 8;
lz     = 2;
numx2  = 6;      % 2^numx2 = number of elements along X direction
numy2  = 6;      % 2^numy2 = number of elements along Y direction
numz2  = 1;      % 2^numz2 = number of elements along X direction
[mesh] = buildGrid3D(lx,ly,lz,numx2,numy2,numz2);
element= mesh.element;
node   = mesh.node;

% move the beam relative to the background grid
node1(:,2) = node1(:,2) + 3.5;
node1(:,1) = node1(:,1) + mesh.deltax;
node1(:,3) = node1(:,3) + 0.5;

markedYCoord = node1(markedNode,2);

% store the particle mesh into a structure for convenience
particles.node     = node1;
particles.elem     = element1;
particles.elemType = elemType;

% find boundary nodes

fixNodes=find(abs(node(:,1)-mesh.deltax)<1e-10);

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum = zeros(mesh.nodeCount,3);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,3);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,3);  % nodal external force vector
nforce    = zeros(mesh.nodeCount,3);  % nodal force vector

%% plot mesh, particles

hold on
plot_mesh(particles.node,particles.elem,elemType,'black-',1.2);
plot_mesh(node,element,'B8','r-',1.2);
plot3(coords(:,1),coords(:,2),coords(:,3),'k.','markersize',10);
plot3(node(fixNodes,1),node(fixNodes,2),node(fixNodes,3),'r*','markersize',14);
axis off

ta = 0;           % time
ka = 0;

% NOTE: sum(mass) = sum(nmass) => particle mass = grid mass

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.4*(mesh.deltax/c);
time  = 3;
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
        data  = getCPDITet4(p,particles,mesh);
        esctr = data.node;
        Vp    = volume(p);
        Mp    = mass(p);
        for i=1:length(esctr)
            id              = esctr(i);
            nmass(id)       = nmass(id)       + data.phi(i)*Mp;
            nmomentum(id,:) = nmomentum(id,:) + data.phi(i)*Mp*velo(p,:);
               
            dNIdx           = data.dphi(i,1);
            dNIdy           = data.dphi(i,2);
            dNIdz           = data.dphi(i,3);
            
            niforce(id,1)   = niforce(id,1) - Vp*(sig(1)*dNIdx + sig(6)*dNIdy + sig(5)*dNIdz);
            niforce(id,2)   = niforce(id,2) - Vp*(sig(6)*dNIdx + sig(2)*dNIdy + sig(4)*dNIdz);
            niforce(id,3)   = niforce(id,3) - Vp*(sig(5)*dNIdx + sig(4)*dNIdy + sig(3)*dNIdz);
                
            neforce(id,:)   = neforce(id,:) +Mp*data.phi(i)*bodyf;           
        end
    end
        
    % update nodal momenta   
    nforce    = niforce + neforce;
    nforce   (fixNodes,:)  = 0;
    nmomentum(fixNodes,:)  = 0;
    
    nmomentum = nmomentum + nforce*dtime;
        
    % update particle velocity and position and stresses
    
    % loop over particles
    for p=1:pCount
        Lp    = zeros(3,3);        
        data  = getCPDITet4(p,particles,mesh);
        esctr = data.node;
        for i=1:length(esctr)
            id = esctr(i);
            vI = [0 0 0];
            if nmass(id) > tol
                velo(p,:)  = velo(p,:) + dtime * data.phi(i)*nforce(id,:)/nmass(id);                
                vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
            end
            Lp = Lp + vI'*data.dphi(i,:);         % particle gradient velocity
        end
        
        F          = (I + Lp*dtime)*reshape(deform(p,:),3,3);
        deform(p,:)= reshape(F,1,9);
        volume(p)  = det(F)*volume0(p);
        J          = det(F);
        b          = F*F';
        sigma      = 1/J*( mu*(b-I) + lambda*log(J)*I );
        stress(p,:)  = [sigma(1,1) sigma(2,2)  sigma(3,3) sigma(2,3) sigma(1,3) sigma(1,2)];
    end
    
    for c=1:size(particles.node,1)        
        xc    = particles.node(c,:);
        ec    = point2ElemIndex3D(xc,mesh);
        esctr = element(ec,:);
        for i=1:length(esctr)
            id      = esctr(i);
            x       = xc - node(id,:);
            [N,dNdx]= getMPM3D(x,mesh.deltax,mesh.deltay,mesh.deltaz);
            if nmass(id) > tol                
                xc = xc + dtime*N*nmomentum(id,:)/nmass(id);                
            end            
        end
        particles.node(c,:) = xc;
    end
    
  
    % VTK output
    
    if (  mod(istep,interval) == 0 )
        vtkFile = sprintf('../../results/cpdi/%s%d',vtkFileName,istep);
        data.stress  = [stress zeros(pCount,1)];
        data.pstrain = [];
        data.color   = color;
        data.velo    = velo;
        VTKParticlesCPDI(particles,vtkFile,data);
    end
    
    
    % advance to the next time step
    
    t     = t + dtime;
    istep = istep + 1;
    
    % compute vertical displacement of marked node and store it   
    u  = particles.node(markedNode,2)-markedYCoord;
    ta = [ta;t];
    ka = [ka;u];        
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

% background grid
Ux= zeros(size(node,1),1);
Uy= zeros(size(node,1),1);
sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

VTKPostProcess3d(node,element,'B8',vtkFileName1,...
             [sigmaXX sigmaYY sigmaXY sigmaXX sigmaXX sigmaXX],[Ux Uy Ux]);

%%



%ss=load('femTet4Cantilever82.mat');
ss1=load('cpdiTet4Cantilever82.mat');
ss2=load('cpdiTet4Cantilever178.mat');
ss3=load('cpdiTet4Cantilever491.mat');
ss4=load('cpdiTet4Cantilever.mat');
sss=load('femTet4Cantilever1000.mat'); % smallest element size=0.0534

ss1g=load('cpdiTet4Cantilever82-bg3.mat');
ss2g=load('cpdiTet4Cantilever178-bg3.mat');
ss3g=load('cpdiTet4Cantilever491-bg3.mat');

figure
set(gca,'FontSize',14)
hold on
%plot(ta(1:end),ka(1:end),'b-','LineWidth',2.);
plot(ss1.ta(1:end),ss1.ka(1:end),'cy-','LineWidth',2.2);
plot(ss2.ta(1:end),ss2.ka(1:end),'b-','LineWidth',2.2);
plot(ss3.ta(1:end),ss3.ka(1:end),'r-','LineWidth',2.2);
plot(ss4.ta(1:end),ss4.ka(1:end),'black-','LineWidth',2.2);
plot(sss.ta(1:end),sss.pDisp(1:end),'r--','LineWidth',2.5);

xlabel('Time [s]')
ylabel('Displacement [m]')
legend('CPDI-Tet4, 82 particles','CPDI-Tet4, 178 particles',...
       'CPDI-Tet4, 491 particles','CPDI-Tet4, 1000 particles', ...
       'FEM-Tet4, 1000 elements')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
grid on
axis([0. 3. 0 -3.5])

%%
figure
set(gca,'FontSize',14)
hold on
%plot(ta(1:end),ka(1:end),'b-','LineWidth',2.);
plot(ss1g.ta(1:end),ss1g.ka(1:end),'b-','LineWidth',2.2);
plot(ss2g.ta(1:end),ss2g.ka(1:end),'cys','LineWidth',2.2);
plot(ss3g.ta(1:end),ss3g.ka(1:end),'black-','LineWidth',2.2);
plot(sss.ta(1:end),sss.pDisp(1:end),'r--','LineWidth',2.5);

xlabel('Time [s]')
ylabel('Displacement [m]')
legend('CPDI-Tet4, 82 particles','CPDI-Tet4, 178 particles',...
       'CPDI-Tet4, 491 particles', ...
       'FEM-Tet4, 1000 elements')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
grid on
axis([0. 3. 0 -3.5])

save('cpdiTet4Cantilever491-bg3.mat', 'ta', 'ka');


disp([num2str(toc),'   DONE '])
