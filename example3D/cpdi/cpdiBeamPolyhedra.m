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
% Compliant cantilever beam with extreme deformation.
%
% Vinh Phu Nguyen
% Monash University
% 18 April, 2016 

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
interval     = 20;
vtkFileName  = 'cantileverPolyhedronReg';
vtkFileName1 = '../../results/cpdi/cantileverPolyhedronGrid';

tic;

disp([num2str(toc),'   INITIALISATION '])

%%   particle distribution from a mesh
%

fprintf(1,'....Load the mesh\n');
fileName = 'beam500-centroid.ply'; % Centroidal Voronoi
fileName = 'beam491-reg.ply';    % Regular Voronoi 

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

markedNode1 = find(abs(node(:,3)-4)<1e-10);
markedNode2 = find(abs(node(:,1)-1)<1e-10);
markedNode3 = find(abs(node(:,2)-0)<1e-10);

markedNode  = intersect(intersect(markedNode1,markedNode2),markedNode3);

pCount  = numelem;                       % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
coords  = zeros(pCount,3); % particle coordinates (centroids)
deform  = repmat([1 0 0 0 1 0 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,6);                % stress
strain  = zeros(pCount,6);                % strain
velo    = zeros(pCount,3);
color   = ones (pCount,1);

%% Computational grid

ghostCell=0;
lx     = 2;
ly     = 8;
lz     = 8;
numx2  = 1;      % 2^numx2 = number of elements along X direction
numy2  = 4;      % 2^numy2 = number of elements along Y direction
numz2  = 4;      % 2^numz2 = number of elements along X direction
[mesh] = buildGrid3D(lx,ly,lz,numx2,numy2,numz2);


% move the beam relative to the background grid
node(:,2) = node(:,2) + 3.5;
node(:,3) = node(:,3) + mesh.deltaz;
node(:,1) = node(:,1) + 0.5;

markedYCoord = node(markedNode,2);

% store the particle mesh into a structure for convenience
particles.node     = node;
particles.element  = element; % element connectivity
particles.elem     = elem;    % faces of every element
particles.face     = Elements.face;
%particles.elemType = elemType;

% reuse CPDI functions to save time

pnodes = cell(pCount,1);
pfuncW = cell(pCount,1);
pgradW = cell(pCount,3);

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
res.wf(:)=0;res.wg(:)=0;

% find boundary nodes

fixNodes=find(abs(mesh.node(:,3)-mesh.deltaz)<1e-6);

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
plot3(coords(:,1),coords(:,2),coords(:,3),'k.','markersize',10);
plot3(mesh.node(fixNodes,1),mesh.node(fixNodes,2),mesh.node(fixNodes,3),'b*','markersize',14);
axis off

ta = 0;           % time
ka = 0;

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.2*(mesh.deltaz/c);
time  = 2.6;
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
        
%         if ( abs(sum(data.phi)-1) > 1e-3 )
%             error('PoU not satisfied!!!') 
%         end;
        
        
        esctr = data.node;
        Vp    = volume(p);
        Mp    = mass(p);
        vp    = velo(p,:);
        
        pnodes{p} = esctr;
        pfuncW{p} = data.phi;
        pgradW{p} = data.dphi;
        
        for i=1:length(esctr)
            id              = esctr(i);
            Ni              = data.phi(i);
            nmass(id)       = nmass(id)       + Ni*Mp;
            nmomentum(id,:) = nmomentum(id,:) + Ni*Mp*vp;
                           
            dNIdx           = data.dphi(i,1);
            dNIdy           = data.dphi(i,2);
            dNIdz           = data.dphi(i,3);
            
            niforce(id,1)   = niforce(id,1) - Vp*(sig(1)*dNIdx + sig(6)*dNIdy + sig(5)*dNIdz);
            niforce(id,2)   = niforce(id,2) - Vp*(sig(6)*dNIdx + sig(2)*dNIdy + sig(4)*dNIdz);
            niforce(id,3)   = niforce(id,3) - Vp*(sig(5)*dNIdx + sig(4)*dNIdy + sig(3)*dNIdz);
                
            neforce(id,:)   = neforce(id,:) + Mp*Ni*bodyf;           
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
        %data  = getCPDIPolyhedronBasis(p,coords(p,:),particles,mesh);
        %esctr = data.node;
        esctr     = pnodes{p};
        data.phi  = pfuncW{p};
        data.dphi = pgradW{p};
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
    
    % Should we have to update centroids???
    for e = 1:numelem
        corners     = particles.node(element{e},:);
        %faces       = elem{e}.faces;
        %coords(e,:) = centroid(corners); % center of each element=particle
        coords(e,:) = mean(corners); % center of each element=particle
        %[res,V]     = getCPDIPolyhedronData(faces,coords(e,:),particles,res);        
        %volume(e)   = V;
    end

    % VTK output
    
    if (  mod(istep,interval) == 0 )
        vtkFile = sprintf('../../results/cpdi/%s%d',vtkFileName,istep);        
        data.stress  = [stress zeros(pCount,1)];        
        data.velo    = velo;
        data.color   = color;
         exportToParaviewPoly(particles.node,glbFace,3,sumc,totalf,...
          zeros(numnode,1),zeros(numnode,2),zeros(numnode,3),vtkFile);
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
Ux= zeros(size(mesh.node,1),1);
Uy= zeros(size(mesh.node,1),1);
sigmaXX = zeros(size(mesh.node,1),1);
sigmaYY = zeros(size(mesh.node,1),1);
sigmaXY = zeros(size(mesh.node,1),1);

VTKPostProcess3d(mesh.node,mesh.element,'B8',vtkFileName1,...
             [sigmaXX sigmaYY sigmaXY sigmaXX sigmaXX sigmaXX],[Ux Uy Ux]);
         
         
% pvdFile = fopen(strcat('../../results/cpdi/','cantileverPolyhedra.pvd'), 'wt');
% 
% fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
% fprintf(pvdFile,'<Collection>\n');
% 
% for i = 1:nsteps
%     if (  mod(i,interval) == 0 )
%         vtuFile = sprintf('%s%d%s',vtkFileName,i,'.ply');
%         fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
%     end
% end
% 
% fprintf(pvdFile,'</Collection>\n');
% fprintf(pvdFile,'</VTKFile>\n');
% 
% fclose(pvdFile);         

%%



ss1 = load('cpdiPolyhedralCantilever491.mat');
ss2 = load('cpdiPolyhedralCantilever300.mat');
%ss3 = load('cpdiTet4Cantilever491.mat');
ss4 = load('femTet4Cantilever1000.mat');

figure
set(gca,'FontSize',14)
hold on
plot(ss1.ta(1:end),ss1.ka(1:end),'black-','LineWidth',2);
plot(ta(1:end),ka(1:end),'b-','LineWidth',2);
%plot(ss2.ta(1:end),ss2.ka(1:end),'black-','LineWidth',2);
%plot(ss3.ta(1:end),ss3.ka(1:end),'cy-','LineWidth',2.3);
plot(ss4.ta(1:end),ss4.pDisp(1:end),'r-','LineWidth',2.3);
xlabel('Time [s]')
ylabel('Displacement [m]')
legend('CPDI-Polyhedron 491, CVT','CPDI-polyhedron 491', 'FEM')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%axis([0 0.25 0 -1800])
grid on
axis([0. 3. 0 -3.5])

save('cpdiPolyhedralCantilever491-reg.mat');

disp([num2str(toc),'   DONE '])
