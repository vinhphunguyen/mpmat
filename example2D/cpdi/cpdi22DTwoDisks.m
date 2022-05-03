% This file implements the Material Point Method with CPDI2 interpolation
% described in the article
%
% A. Sadeghirad, R. M. Brannon, and J.E. Guilkey. Second-order convected
% particle domain in- terpolation (CPDI2) with enrichment for weak discontinuities
% at material interfaces. IJNME, 95(11):928-952, 2013.
%
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% Particles taken from a Q4 unstructured mesh (Gmsh).
%
% Two elastic disks come into contact.
% Example taken from Sulsky et al. 1994 paper.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% May 2014, Saigon, Vietnam.

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
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%

E   = 1000;        % Young's modulus
nu  = 0.3;         % Poisson ratio
rho = 1000;        % density
kappa = 3-4*nu;    % Kolosov constant
mu    = E/2/(1+nu);% shear modulus

v   = 0.1;

interval     = 100;
vtkFileName  = 'cpdi2DTwoDisks';
vtkFileName1 = '../results/cpdi2/cpdi2DTwoDisksGrid';

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell=1;
l     = 1;
numx2 = 20;      % number of elements along X direction
numy2 = 20;      % number of elements along Y direction
[mesh]= buildGrid2D(l,l,numx2,numx2, ghostCell);
element= mesh.element;
node   = mesh.node;

%%   particle distribution from a mesh
%

meshFile = 'disks_q4.msh';
pmesh     = load_gmsh (meshFile);

elemType = 'Q4';
numnode  = pmesh.nbNod;
numelem  = pmesh.nbQuads;
node1    = pmesh.POS(:,1:2);
element1 = pmesh.QUADS(1:numelem,1:4);

% move the particles due to ghostCell
node1(:,1) = node1(:,1) + ghostCell*mesh.deltax;
node1(:,2) = node1(:,2) + ghostCell*mesh.deltay;

% store the particle mesh into a structure for convenience
particles.node     = node1;
particles.elem     = element1;
particles.elemType = elemType;

pCount  = numelem;                       % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
coords  = zeros(pCount,2);
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stressx
strain  = zeros(pCount,3);                % strain
velo    = zeros(pCount,2);
nodeid  = cell(pCount,1);
funcW   = zeros(pCount,4);
gradW   = zeros(pCount,8);

for e = 1:numelem
    coord = node1(element1(e,:),:);
    a     = 0.5*( coord(1,1)*coord(2,2)  - coord(2,1)*coord(1,2) ...
                + coord(2,1)*coord(3,2)  - coord(3,1)*coord(2,2) ...
                + coord(3,1)*coord(4,2)  - coord(4,1)*coord(3,2) ...
                + coord(4,1)*coord(1,2)  - coord(1,1)*coord(4,2) );
    volume(e)  = a;
    mass(e)    = a*rho;
    coords(e,:) = mean(coord); % center of each element=particle
    if coords(e,1) < 0.5
        velo(e,:) = [v v];
    else
        velo(e,:) = [-v -v];
    end
end

volume0 = volume;

for p=1:pCount
  data = getCPDIQuadData(p,particles,mesh);
  nodeid{p}  = data.nodes;
  funcW(p,:) = data.wf;
  gradW(p,:) = data.wg(:);
end

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum = zeros(mesh.nodeCount,2);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,2);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,2);  % nodal external force vector

%% plot mesh, particles

hold on
plot_mesh(particles.node,particles.elem,elemType,'b-',2.1);
%plot_mesh(particles.node,particles.elem(1,:),elemType,'r--',0.2);
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(coords(:,1),coords(:,2),'k.','markersize',30);
axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.001;
time  = 3.5;
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
    disp(['time step ',num2str(t)])
    % reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    % loop over particles
    for p=1:pCount
        sig    = stress(p,:);
        % particle mass and momentum to node
        %data  = cpdi22D(p,particles,mesh);
        input.nodes=nodeid{p};
        input.wf   = funcW(p,:);
        input.wg   = reshape(gradW(p,:),4,2);
        input.Vp   = volume(p);
        
        data  = getCPDIQuadBasis(p,input,particles,mesh);
        esctr = data.node;
        for i=1:length(esctr)
            id              = esctr(i);
            nmass(id)       = nmass(id)       + data.phi(i)*mass(p);
            nmomentum(id,:) = nmomentum(id,:) + data.phi(i)*mass(p)*velo(p,:);
            niforce(id,1)   = niforce(id,1) - volume(p)*(sig(1)*data.dphi(i,1) + sig(3)*data.dphi(i,2));
            niforce(id,2)   = niforce(id,2) - volume(p)*(sig(3)*data.dphi(i,1) + sig(2)*data.dphi(i,2));
        end
    end
    
    % debug
    % update nodal momenta
    
    nmomentum = nmomentum + niforce*dtime;
    % no bounda
    % update particle velocity and position and stresses
    k = 0;
    u = 0;
    % loop over particles
    for p=1:pCount
        Lp    = zeros(2,2);        
        %data  = cpdi22D(p,particles,mesh); % old implementation
        input.nodes=nodeid{p};
        input.wf   = funcW(p,:);
        input.wg   = reshape(gradW(p,:),4,2);
        input.Vp   = volume(p);
        
        data  = getCPDIQuadBasis(p,input,particles,mesh);
        esctr = data.node;
        for i=1:length(esctr)
            id = esctr(i);
            vI = [0 0];
            if nmass(id) > tol
                velo(p,:)  = velo(p,:) + dtime * data.phi(i)*niforce(id,:)  /nmass(id);
                %coords(p,:)= coords(p,:)+ dtime * data.phi(i)*nmomentum(id,:)/nmass(id);
                vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
            end
            Lp = Lp + vI'*data.dphi(i,:);         % particle gradient velocity
        end
        
        F          = ([1 0;0 1] + Lp*dtime)*reshape(deform(p,:),2,2);
        deform(p,:)= reshape(F,1,4);
        volume(p)  = det(F)*volume0(p);
        dEps       = dtime * 0.5 * (Lp+Lp');
        dsigma     = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
        stress(p,:)= stress(p,:) + dsigma';
        strain(p,:)= strain(p,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
        
        k = k + 0.5*(velo(p,1)^2+velo(p,2)^2)*mass(p);
        %             u = u + 0.25/mu*( 0.25*(kappa+1)*(s(pid,1)^2+s(pid,2)^2) ...
        %                   - 2*(s(pid,1)*s(pid,2)-s(pid,3)^2))*Vp(pid);
        u = u + 0.5*volume(p)*stress(p,:)*strain(p,:)';
    end
    
    % update particle corners position
    for c=1:numnode        
        xc    = particles.node(c,:);
        ec    = point2ElemIndex(xc,mesh);
        esctr = element(ec,:);
        for i=1:length(esctr)
            id = esctr(i);
            x     = xc - node(id,:);
            [N,dNdx]=getMPM2D(x,mesh.deltax,mesh.deltay);
            if nmass(id) > tol                
                xc = xc + dtime*N*nmomentum(id,:)/nmass(id);                
            end            
        end
        particles.node(c,:) = xc;
    end
    
    % update CPDI data 
    
    for p=1:pCount
      data = getCPDIQuadData(p,particles,mesh);
      nodeid{p}  = data.nodes;
      funcW(p,:) = data.wf;
      gradW(p,:) = data.wg(:);
    end
    
    % store time,velocty for plotting
    
    ta = [ta;t];
    ka = [ka;k];
    sa = [sa;u];
    
    % VTK output
    
    if (  mod(istep,interval) == 0 )
        vtkFile = sprintf('../results/cpdi2/%s%d',vtkFileName,istep);
        data.stress  = [stress zeros(pCount,1)];
        data.pstrain = [];
        data.velo    = velo;
        VTKParticlesCPDI(particles,vtkFile,data);
    end
    
    
    % advance to the next time step
    
    t     = t + dtime;
    istep = istep + 1;
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

pvdFile = fopen(strcat('../results/cpdi2/',vtkFileName,'.pvd'), 'wt');

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

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(ta(1:end),sa(1:end),'r--','LineWidth',2);
plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy')
legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
axis([0 3 0 3])

ss=load('gimp-2disks-energies.mat');

figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'g-','LineWidth',2.1);
plot(ss.ta(1:end),ss.ka(1:end),'r-','LineWidth',2.1);
plot(ta(1:end),sa(1:end),'g--','LineWidth',2.1);
plot(ss.ta(1:end),ss.sa(1:end),'r--','LineWidth',2.1);
xlabel('Time')
ylabel('Energy')
legend('cpdi','gimp')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%%
figure
hold on
plot_mesh(particles.node,particles.elem,elemType,'r-',2.);
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(coords(:,1),coords(:,2),'k.','markersize',10);
axis off

disp([num2str(toc),'   DONE '])
