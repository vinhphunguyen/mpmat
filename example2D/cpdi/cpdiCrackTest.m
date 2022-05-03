% This file implements the Material Point Method with CPDI-Q4 interpolation
% described in the article
%
% A. Sadeghirad, R. M. Brannon, and J.E. Guilkey. Second-order convected
% particle domain in- terpolation (CPDI2) with enrichment for weak discontinuities
% at material interfaces. IJNME, 95(11):928-952, 2013.
%
% Simple test to see how cracks can be modeled with CPDI-Q4.
% Idea: CPDI mesh with duplicated nodes to allow each particle is
% independent from others.
%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 19 August 2015.

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

E     = 1000;        % Young's modulus
nu    = 0.0;         % Poisson ratio: make it 1D
rho   = 1000;        % density
kappa = 3-4*nu;      % Kolosov constant
mu    = E/2/(1+nu);  % shear modulus
ft    = 100.0;       % tensile strength

v0    = 0.2;         % imposed velocity 

interval     = 1;
vtkFileName  = 'cpdiCrack';
vtkFileName1 = '../results/cpdi2/crack/cpdiCrackGrid';

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C           = elasticityMatrix(E,nu,stressState);
D           = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell=1;
l     = 7;
numx2 = 7;      % number of elements along X direction
numy2 = 7;      % number of elements along Y direction
[mesh]= buildGrid2D(l,l,numx2,numx2, ghostCell);
mesh.node = mesh.node - [ones(mesh.nodeCount,1)*mesh.deltax ones(mesh.nodeCount,1)*mesh.deltay];
element= mesh.element;
node   = mesh.node;

% find boundary nodes

eps=1e-12;
leftNodes  = mesh.lNodes;
rightNodes = mesh.rNodes;

%%   particle distribution from a mesh
%

node1 = [1.5 2.5;
         5.5 2.5;
         5.5 4.5;
         1.5 4.5;
         3.5 2.5; 
         3.5 4.5; 
         3.5 2.5; % extra nodes for cracking
         3.5 4.5; % extra nodes for cracking
         1.01   2.5;
         1.01   4.5;
         5.99   2.5;
         5.99   4.5
        ];

     
element1    = cell(4,1);
element1{1} = [1 5 6 4]; 
element1{2} = [7 2 3 8];      
element1{3} = [9 1 4 10];
element1{4} = [2 11 12 3];
          
numelem  = size(element1,1);
numnode  = size(node1,1); 
elemType = 'Q4';

% store the particle mesh into a structure for convenience
particles.node     = node1;
particles.elem     = element1;
particles.elemType = elemType;

pCount  = numelem;                       % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
position= zeros(pCount,2);                % particle position
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stressx
strain  = zeros(pCount,3);                % strain
velo    = zeros(pCount,2);                % velo
color   = zeros(pCount,1);                % color
crack   = zeros(pCount,1);                % cracked or not
% CPDI data
nodeid  = cell(pCount,1);                 % nodes affect particle 'p'
funcW   = cell(pCount,1);                 % function weights of 'p'
gradW   = cell(pCount,2);                 % gradient weights of 'p', 
volum   = cell(pCount,1);                 % store particle domain volume

% color to differentiate particles in post-processing
color(1) = 1;
color(2) = 2;
color([3 4]) = 3;


%crack(:) = 1;

% particle mass, volume and initial velocity
for e = 1:numelem
    coord = node1(element1{e},:);
    a     = 0.5*( coord(1,1)*coord(2,2)  - coord(2,1)*coord(1,2) ...
                + coord(2,1)*coord(3,2)  - coord(3,1)*coord(2,2) ...
                + coord(3,1)*coord(4,2)  - coord(4,1)*coord(3,2) ...
                + coord(4,1)*coord(1,2)  - coord(1,1)*coord(4,2) );
    volume(e)   = a;
    mass(e)     = a*rho;
    position(e,:) = mean(coord); % center of each element=particle
end

volum{1} = volume(1);
volum{2} = volume(2);

volume0 = volume;

for p=1:pCount
  data       = getCPDIQuadDataGeneral(p,particles,mesh);
  nodeid{p}  = data.nodes;
  funcW{p}   = data.wf;
  gradW{p,1} = data.wg(:,1);
  gradW{p,2} = data.wg(:,2);
end

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum0= zeros(mesh.nodeCount,2);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,2);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,2);  % nodal external force vector

%% plot mesh, particles

figure
hold on
plot_mesh(particles.node,cell2mat(particles.elem),elemType,'r-',0.2);
plot_mesh(node,element,'Q4','k-',1.); % background grid
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(position(:,1),position(:,2),'k.','markersize',10);
plot(node(leftNodes,1),node(leftNodes,2),'b*','markersize',10);
plot(node(nodeid{1},1),node(nodeid{1},2),'rs','markersize',10);
axis on % axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.09*mesh.deltax/c;
time  = 1.5;
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
    disp(['time step ',num2str(t)])
    % reset grid data
    nmass(:)     = 0;
    nmomentum0(:)= 0;
    niforce(:)   = 0;
    % loop over deformable particles
    for p=1:2
        sig    = stress(p,:);
        % particle mass and momentum to node
        %data  = cpdi22D(p,particles,mesh);
        input.nodes=nodeid{p};
        input.wf   = funcW{p};
        input.wg   = [gradW{p,1} gradW{p,2}];
        input.Vp   = volume(p);
        
        data  = getCPDIQuadBasisGeneral(p,input,particles,mesh);
        esctr = data.node;
        for i=1:length(esctr)
            id              = esctr(i);
            nmass(id)       = nmass(id)       + data.phi(i)*mass(p);
            nmomentum0(id,:) = nmomentum0(id,:) + data.phi(i)*mass(p)*velo(p,:);
            niforce(id,1)   = niforce(id,1) - volume(p)*(sig(1)*data.dphi(i,1) + sig(3)*data.dphi(i,2));
            niforce(id,2)   = niforce(id,2) - volume(p)*(sig(3)*data.dphi(i,1) + sig(2)*data.dphi(i,2));
        end
    end
    
    % debug
    % update nodal momenta
    
    % boundary conditions on left/right edges of the box
%     nmomentum(leftNodes ,1) = 0;
%     nmomentum(rightNodes,1) = 0;
%     niforce  (leftNodes ,1) = 0;
%     niforce  (rightNodes,1) = 0;
    
    nmomentum = nmomentum0 + niforce*dtime;
    
    
    % boundary conditions from rigid particles
    
    nmomentum(nodeid{3},2) = 0; 
    nmomentum(nodeid{4},2) = 0;
    nmomentum(nodeid{3},1) = nmass(nodeid{3})*(-v0);
    nmomentum(nodeid{4},1) = nmass(nodeid{4})*( v0);
%     niforce(nodeid{3},:)   = 0; 
%     niforce(nodeid{4},:)   = 0; 
    
    % update particle velocity and position and stresses
    k = 0;
    u = 0;
    % loop over deformable particles
    for p=1:2
        Lp    = zeros(2,2);        
        %data  = cpdi22D(p,particles,mesh); % old implementation
        input.nodes=nodeid{p};
        input.wf   = funcW{p};
        input.wg   = [gradW{p,1} gradW{p,2}];
        input.Vp   = volume(p);
        
        data  = getCPDIQuadBasisGeneral(p,input,particles,mesh);
        esctr = data.node;
        for i=1:length(esctr)
            id = esctr(i);
            vI = [0 0];
            if nmass(id) > tol
                velo(p,:)  = velo(p,:) + data.phi(i)*(nmomentum(id,:) - nmomentum0(id,:)) /nmass(id);
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
        
        % simple erosion technique
        %if crack(p), stress(p,:) = 0; end
        
        k = k + 0.5*(velo(p,1)^2+velo(p,2)^2)*mass(p);
        u = u + 0.5*volume(p)*stress(p,:)*strain(p,:)';
        
        % check failure (not sure where)
%         
        if stress(p,1) > ft 
          disp('cracking ..')
          crack(p) = 1;
        end
    end
    
    % update particle corners position
    for p=1:pCount                          % loop over particles
      pCorners = particles.elem{p};         % corners of particle 'p'
      if crack(p) == 0
        [particles] = updateCorners (pCorners,particles,mesh,...
                                     nmomentum,nmass,dtime);
         position(p,:) = mean ( particles.node(pCorners,:) );                                   
      else
        enode = particles.node(pCorners,:);
        pt(1) = (2*position(p,1)-(enode(1,1)+enode(2,1)))/(enode(2,1)-enode(1,1));
        pt(2) = (2*position(p,2)-(enode(2,2)+enode(3,2)))/(enode(3,2)-enode(2,2));
        
        [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
      
        for i=1:length(N)
          cid = pCorners(i);    
          xc  = particles.node(cid,:) + dtime*N(i)*velo(p,:);        
          particles.node(cid,:) = xc;
        end  
        position(p,:) = mean ( particles.node(pCorners,:) );                                   
      end      
    end
    
    
    % update CPDI data 
    
    for p=1:pCount
      data = getCPDIQuadDataGeneral(p,particles,mesh);
      nodeid{p}  = data.nodes;
      funcW{p}   = data.wf;
      gradW{p,1} = data.wg(:,1);
      gradW{p,2} = data.wg(:,2);
      volum{p}   = data.Vp; 
      volume(p)=sum(data.Vp);
    end
    
    % store time,velocty for plotting
    
    ta = [ta;t];
    ka = [ka;k];
    sa = [sa;u];
    
    % VTK output
    
    if (  mod(istep,interval) == 0 )
        vtkFile = sprintf('../results/cpdi2/crack/%s%d',vtkFileName,istep);
        data.stress  = [stress zeros(pCount,1)];
        data.pstrain = [];
        data.velo    = velo;
        data.color   = color;
        VTKParticlesCPDI(particles,vtkFile,data);
    end
    
    
    % advance to the next time step
    
    t     = t + dtime;
    istep = istep + 1;
end



%% post processing

% disp([num2str(toc),'   POST-PROCESSING '])
% 
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
Ux= zeros(size(mesh.node,1),1);
Uy= zeros(size(mesh.node,1),1);
sigmaXX = zeros(size(mesh.node,1),1);
sigmaYY = zeros(size(mesh.node,1),1);
sigmaXY = zeros(size(mesh.node,1),1);

VTKPostProcess(mesh.node,mesh.element,2,'Quad4',vtkFileName1,...
    [sigmaXX sigmaYY sigmaXY],[Ux Uy]);

%%

% ss=load('cpdi-hole-normal-energies.mat');
% 
figure
subplot(2,1,1)
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'black-','LineWidth',2.1);
plot(ta(1:end),sa(1:end),'red-','LineWidth',2.1);
xlabel('Time')
ylabel('Kinetic energy')



figure
hold on
plot_mesh(particles.node,cell2mat(particles.elem),elemType,'r-',0.2);
plot_mesh(node,element,'Q4','k-',1.); % background grid
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(position(:,1),position(:,2),'k.','markersize',10);
axis off

disp([num2str(toc),'   DONE '])
