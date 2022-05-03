% This file implements the Material Point Method with CPDI-Q4 interpolation
% described in the article
%
% Demonstration of polygonal CPDI in MPM.
%
% Polygonal mesh is generated using PolyMesher developed by Paulino et al.
% Using the sub-sampling method to derive n-gon CPDI functions.
%
%
% Vinh Phu Nguyen
% The University of Adelaide, SA, Australia
% 22 August 2015.

%%

addpath ../../fem_util/
addpath ../../fem-functions/
addpath ../../post-processing/
addpath ../../externals/PolyMesher/

%%
clc
clear 
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

v0    = 2;         % imposed velocity 

interval     = 1;
vtkFileName  = 'cpdiPolygon';
vtkFileName1 = '../results/cpdi2/polygon/cpdiGrid';

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C           = elasticityMatrix(E,nu,stressState);
D           = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

ghostCell=1;
l     = 0.7;
numx2 = 25;      % number of elements along X direction
numy2 = 25;      % number of elements along Y direction
[mesh]= buildGrid2D(l,l,numx2,numx2, ghostCell);
mesh.node = mesh.node - [ones(mesh.nodeCount,1)*mesh.deltax ones(mesh.nodeCount,1)*mesh.deltay];
element= mesh.element;
node   = mesh.node;

% find boundary nodes

eps=1e-12;
leftNodes  = mesh.lNodes;
botNodes   = mesh.bNodes;

%%   particle distribution from a mesh
% generate Voronoi mesh where CircleDomain.m defines the geometry
% which is in this case a circle
NElem = 20;
[Node,Element,Supp,Load,P]=PolyMesher(@CircleDomain,NElem,100);

% store the particle mesh into a structure for convenience
particles.node     = Node;
particles.elem     = Element;

pCount  = size(Element,1);                % # of particles
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

% initial particle velocity

velo(:,2) = -v0;

% particle mass, volume and initial velocity
for e = 1:NElem
  vx  = Node(Element{e},1); 
  vy  = Node(Element{e},2); 
  nv  = length(Element{e});
  vxS = vx([2:nv 1]); 
  vyS = vy([2:nv 1]); %Shifted vertices
  temp        = vx.*vyS - vy.*vxS;
  a           = 0.5*sum(temp);
  volume(e)   = a;
  mass(e)     = a*rho;
  position(e,:) = 1/(6*a)*[sum((vx+vxS).*temp),sum((vy+vyS).*temp)];
end

volume0 = volume;

for p=1:pCount
  data       = getCPDIPolygonData(p,particles,mesh);
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
MaxNVer = max(cellfun(@numel,Element));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Element,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Node,'FaceColor','w'); pause(1e-6)
%plot_mesh(node,element,'Q4','k-',1.); % background grid
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(position(:,1),position(:,2),'k.','markersize',10);
%plot(node(leftNodes,1),node(leftNodes,2),'b*','markersize',10);
axis on % axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.09*mesh.deltax/c;
time  = 5*dtime;
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < 0 )
    disp(['time step ',num2str(t)])
    % reset grid data
    nmass(:)     = 0;
    nmomentum0(:)= 0;
    niforce(:)   = 0;
    % loop over deformable particles
    for p=1:pCount
      sig    = stress(p,:);
      % particle mass and momentum to node
      input.nodes=nodeid{p};
      input.wf   = funcW{p};
      input.wg   = [gradW{p,1} gradW{p,2}];
      input.Vp   = volume(p);
      
      shape = getCPDIPolygonBasis(p,input,particles,mesh);
      esctr = shape.node;
      for i=1:length(esctr)
        id              = esctr(i);
        nmass(id)       = nmass(id)        + shape.phi(i)*mass(p);
        nmomentum0(id,:)= nmomentum0(id,:) + shape.phi(i)*mass(p)*velo(p,:);
        niforce(id,1)   = niforce(id,1) - volume(p)*(sig(1)*shape.dphi(i,1) + sig(3)*shape.dphi(i,2));
        niforce(id,2)   = niforce(id,2) - volume(p)*(sig(3)*shape.dphi(i,1) + sig(2)*shape.dphi(i,2));
      end
    end
    
    % debug
    % update nodal momenta
    
    % boundary conditions on left/right edges of the box
     nmomentum(botNodes ,:) = 0;
%     nmomentum(rightNodes,1) = 0;
     niforce  (botNodes ,:) = 0;
%     niforce  (rightNodes,1) = 0;
    
    nmomentum = nmomentum0 + niforce*dtime;

    % update particle velocity and position and stresses
    k = 0;
    u = 0;
    % loop over deformable particles
    for p=1:pCount
        Lp    = zeros(2,2);        
        %data  = cpdi22D(p,particles,mesh); % old implementation
        input.nodes=nodeid{p};
        input.wf   = funcW{p};
        input.wg   = [gradW{p,1} gradW{p,2}];
        input.Vp   = volume(p);
        
        shape  = getCPDIPolygonBasis(p,input,particles,mesh);
        esctr  = shape.node;
        for i=1:length(esctr)
          id = esctr(i);
          vI = [0 0];
          if nmass(id) > tol
            velo(p,:)  = velo(p,:) + shape.phi(i)*(nmomentum(id,:) - nmomentum0(id,:)) /nmass(id);
            vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
          end
          Lp = Lp + vI'*shape.dphi(i,:);         % particle gradient velocity
        end
        
        F          = ([1 0;0 1] + Lp*dtime)*reshape(deform(p,:),2,2);
        deform(p,:)= reshape(F,1,4);
        volume(p)  = det(F)*volume0(p);
        dEps       = dtime * 0.5 * (Lp+Lp');
        dsigma     = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
        stress(p,:)= stress(p,:) + dsigma';
        strain(p,:)= strain(p,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
        
        k = k + 0.5*(velo(p,1)^2+velo(p,2)^2)*mass(p);
        u = u + 0.5*volume(p)*stress(p,:)*strain(p,:)';
    end
    
    % update particle corners position
    
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
    
    % update CPDI data 
    
    for p=1:pCount
      data       = getCPDIPolygonData(p,particles,mesh);
      nodeid{p}  = data.nodes;
      funcW{p}   = data.wf;
      gradW{p,1} = data.wg(:,1);
      gradW{p,2} = data.wg(:,2);
      volum{p}   = data.Vp; 
      volume(p)  = sum(data.Vp);
    end
    
    % store time,velocty for plotting
    
    ta = [ta;t];
    ka = [ka;k];
    sa = [sa;u];
    
    % VTK output
    
%     if (  mod(istep,interval) == 0 )
%         vtkFile = sprintf('../results/cpdi2/crack/%s%d',vtkFileName,istep);
%         data.stress  = [stress zeros(pCount,1)];
%         data.pstrain = [];
%         data.velo    = velo;
%         data.color   = color;
%         VTKParticlesCPDI(particles,vtkFile,data);
%     end
%     
    
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
plot_mesh(node,element,'Q4','k-',1.); % background grid
patch('Faces',ElemMat,'Vertices',particles.node,'FaceColor','w'); pause(1e-6)
patch('Faces',ElemMat,'Vertices',Node,'FaceColor','b'); pause(1e-6)

disp([num2str(toc),'   DONE '])
