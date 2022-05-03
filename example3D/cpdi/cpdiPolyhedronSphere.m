% This file implements the Material Point Method with CPDI2 interpolation
% described in the article
%
% A. Sadeghirad, R. M. Brannon, and J.E. Guilkey. Second-order convected
% particle domain in- terpolation (CPDI2) with enrichment for weak discontinuities
% at material interfaces. IJNME, 95(11):928-952, 2013.
%
% Three dimensional problems.
% The grid is a structured mesh consisting of 8-noded trilinear elements (H8).
% CPDI-Polyhedra: polyhedra mesh.
%
% Sphere impacts with rigid walls.
%
% Vinh Phu Nguyen
% Monash University
% 29 April, 2016 

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
v0     = 200;

I = [1 0 0;0 1 0; 0 0 1]; % identity matrix

% be careful with vtkFileName1 and change it according to your computer!!!
interval     = 1;
vtkFileName  = 'spherePolyhedron';
vtkFileName1 = '../../results/cpdi/spherePolyhedronGrid';

tic;

disp([num2str(toc),'   INITIALISATION '])

%%   particle distribution from a mesh
%

fprintf(1,'....Load the mesh\n');
fileName = 'sphere200.ply';

[node,element,Elements] = getPolyMeshdata(fileName);
numelem = size(element,2);
numnode = size(node,1);

% compute total number of cells...
% total number of cells = 
sumc = 0;
for iel = 1:size(element,2)
    nec = element{iel};
    sumc = sumc + (1+length(nec));
end
 
for im = 1:size(Elements.face.vertex_indices,1)
    glbFace{im} = Elements.face.vertex_indices{im} + 1 ;
end
 
% total number of faces...
totalf = 0;
for iel = 1:size(glbFace,2)
    nec = glbFace{iel};
    totalf = totalf + (1 + length(nec));
end

%... vertex indices for all the polytopes...
Faces = Elements.face.vertex_indices;

% loop over elements...
for iel = 1:numelem
    
    % face ids of the current element....
    faceid = Elements.cell.face_indices{iel} + 1;
    
    clear faces
    nds = [];
    % get the faces connecting this particular element
    for in = 1:length(faceid)
        
        % corner indexes of the current polytope
        tm = Faces{faceid(in)} + 1;
        
        % nodal coordinates of the current face
        nds = [nds; node(tm,:)] ;
        
        faces{in} = tm;
    end
    
    %...this saves the list of corner indexes for the current polytope
    elem{iel,1}.faces = faces;
end

% marked corner (recording displacement)

markedNode1 = find(abs(node(:,1)-1)<1e-10);
markedNode2 = find(abs(node(:,2)-0)<1e-10);
markedNode3 = find(abs(node(:,3)-4)<1e-10);

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

velo(:,2) = v0;

%% Computational grid

ghostCell=0;
lx     = 6;
ly     = 6;
lz     = 6;
numx2  = 4;      % 2^numx2 = number of elements along X direction
numy2  = 4;      % 2^numy2 = number of elements along Y direction
numz2  = 4;      % 2^numz2 = number of elements along X direction
[mesh] = buildGrid3D(lx,ly,lz,numx2,numy2,numz2);


% move the beam relative to the background grid
node(:,2) = node(:,2) + 3;
node(:,3) = node(:,3) + 3;
node(:,1) = node(:,1) + 3;


% store the particle mesh into a structure for convenience
particles.node     = node;
particles.element  = element; % element connectivity
particles.elem     = elem;    % faces of every element
%particles.elemType = elemType;


res.wf = zeros(numnode+1,1); res.wg = zeros(numnode+1,3);
for e = 1:numelem
    corners     = node(element{e},:);
    faces       = elem{e}.faces;
    %coords(e,:) = centroid(corners); % center of each element=particle  
    coords(e,:) = mean(corners); % center of each element=particle 
    [res,V]     = getCPDIPolyhedronData(faces,coords(e,:),particles,res);
   
    volume(e)   = V;
    mass(e)     = V*rho;    
end

volume0 = volume;


% find boundary nodes

fixNodes=[find(abs(mesh.node(:,2)-0)<1e-10);...
          find(abs(mesh.node(:,2)-6)<1e-10)];

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum = zeros(mesh.nodeCount,3);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,3);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,3);  % nodal external force vector
nforce    = zeros(mesh.nodeCount,3);  % nodal force vector

%% plot mesh, particles

hold on
%plot_mesh(particles.node,particles.elem,elemType,'black-',1.2);
plot_mesh(mesh.node,mesh.element,'B8','r-',1.2);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
plot(mesh.node(fixNodes,1),mesh.node(fixNodes,2),'r*','markersize',14);
axis off

ta = 0;           % time
ka = 0;

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.05*(mesh.deltax/c);
time  = 20*dtime;
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
        data  = getCPDIPolyhedronBasis(p,coords(p,:),particles,mesh);
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
        data  = getCPDIPolyhedronBasis(p,coords(p,:),particles,mesh);
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
        esctr = mesh.element(ec,:);
        for i=1:length(esctr)
            id      = esctr(i);
            x       = xc - mesh.node(id,:);
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
        exportToParaviewPoly(particles.node,glbFace,3,sumc,totalf,...
          velo(:,1),velo(:,2),velo(:,3),vtkFile);
    end
    
    
    % advance to the next time step
    
    t     = t + dtime;
    istep = istep + 1;
        
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

% background grid
Ux= zeros(size(mesh.node,1),1);
Uy= zeros(size(mesh.node,1),1);
sigmaXX = zeros(size(mesh.node,1),1);
sigmaYY = zeros(size(mesh.node,1),1);
sigmaXY = zeros(size(mesh.node,1),1);

VTKPostProcess3d(mesh.node,mesh.element,'B8',vtkFileName1,...
             [sigmaXX sigmaYY sigmaXY sigmaXX sigmaXX sigmaXX],[Ux Uy Ux]);

%%



%ss=load('cpdi2VerticalBar.mat');

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b*','LineWidth',1.6);
xlabel('Time')
ylabel('Displacement')
legend('CPDI2','CPDI')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%axis([0 0.25 0 -1800])

%save('cpdiTet4Cantilever.mat', 'ta', 'ka');


disp([num2str(toc),'   DONE '])
