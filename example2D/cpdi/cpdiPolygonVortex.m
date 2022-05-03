% This file implements the Material Point Method with CPDI-Q4 interpolation
% described in the article
%
% Demonstration of polygonal CPDI in MPM.
%
% Problem: generalized vortex problem.
%
% Polygonal mesh is generated using PolyMesher developed by Paulino et al.
% Using the sub-sampling method to derive n-gon CPDI functions.
%
%
% Vinh Phu Nguyen
% 22 January 2016.
%


%%
addpath ../../grid/
addpath ../../basis/
addpath ../../particleGen/
addpath ../../constitutiveModels/
addpath ../../util/
addpath ../../geoMesh/
addpath ../../externals/PolyMesher/
addpath ../../postProcessing/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%
E      = 1000;    % Young's modulus
nu     = 0.3;     % Poisson ratio
rho    = 1000;    % density
K      = E/3/(1-2*nu);    % bulk modulus
mu     = E/2/(1+nu);% shear modulus
lambda = K - 2/3*mu;
G      = 1;
Ri     = 0.75;
Ro     = 1.25;

I  = [1 0;0 1];

interval     = 1;
vtkFileName  = 'cpdiPolygonVortex';
vtkFileName1 = '../../results/cpdi/cpdiPolyVortexGrid';

intervalRemesh = 5;

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid
lx        = 3;
ly        = 3;
ghostCell = 0;
numx2     = 40;      % number of elements along X direction
numy2     = 40;      % number of elements along Y direction
[mesh]    = buildGrid2D(lx,ly,numx2,numy2, ghostCell);

mesh.node(:,1) = mesh.node(:,1) - lx/2;
mesh.node(:,2) = mesh.node(:,2) - lx/2;
%mesh.node = mesh.node - [ones(mesh.nodeCount,1)*mesh.deltax ones(mesh.nodeCount,1)*mesh.deltay];

element   = mesh.element;
node      = mesh.node;
elemCount = mesh.elemCount;
nodeCount = mesh.nodeCount;
deltax    = mesh.deltax;
deltay    = mesh.deltay;

% find nodes not in the ring

idx=[];
for i=1:nodeCount
    xI  = node(i,1);
    yI  = node(i,2);
    r   = sqrt(xI*xI+yI*yI);
    if ( ( r > Ro ) || ( r < Ri ) )
        idx = [idx;i];
    end
end

%%   particle distribution from a mesh
% generate Voronoi mesh where RingDomain.m defines the geometry
% which is in this case a circle
NElem = 600;
[Node,Element,Supp,Load,P]=PolyMesher(@RingDomain,NElem,30);

% % write the polygonal elements to file 
% 
% file = 'twoDisks.poly';
% writePolyMesh(file,Node,Element);

% store the particle mesh into a structure for convenience
particles.node     = Node;
particles.elem     = Element;

pCount  = size(Element,1);                % # of particles
volume  = zeros(pCount,1);
mass    = zeros(pCount,1);
position= zeros(pCount,2);                % particle position (centroid of polygons)
deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
stress  = zeros(pCount,3);                % stressx
strain  = zeros(pCount,3);                % strain
velo    = zeros(pCount,2);                % velo
u       = zeros(pCount,2);                % displacement
color   = zeros(pCount,1);                % color
% CPDI data
nodeid  = cell(pCount,1);                 % nodes affect particle 'p'
funcW   = cell(pCount,1);                 % function weights of 'p'
gradW   = cell(pCount,2);                 % gradient weights of 'p', 
volum   = cell(pCount,1);                 % store particle domain volume

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

position0 = position;
volume0   = volume;

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
nvelo     = zeros(mesh.nodeCount,2);  % nodal external force vector

%% New: store shape functions/grads/nodes for all particles
% to speed up 

basis    = cell(pCount,1);
grad     = cell(pCount,1);

%% plot mesh, particles

figure
hold on
MaxNVer = max(cellfun(@numel,Element));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Element,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Node,'FaceColor','w'); 
plot_mesh(node,element,'Q4','k-',1.); % background grid
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(position(end,1),position(end,2),'k.','markersize',10);
plot(node(idx,1),node(idx,2),'b*','markersize',10);

axis on % axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.05*(mesh.deltax/c);
time  = 1; %time=9e-06;
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
    niforce(:)   = 0; neforce(:)   = 0;
    nvelo(:)     = 0;
    % loop over deformable particles
    for p=1:pCount
      sig    = stress(p,:);
      xp     = position(p,:);
      Xp     = position0(p,:);
      mp     = mass(p);
      volp   = volume(p);
      % particle mass and momentum to node
      input.nodes=nodeid{p};
      input.wf   = funcW{p};
      input.wg   = [gradW{p,1} gradW{p,2}];
      input.Vp   = volume(p);
      [bx,by]    = vortexBodyForces(Xp(1),Xp(2),t,mu,rho,G,Ri,Ro);
      shape = getCPDIPolygonBasis(p,input,particles,mesh);
      esctr = shape.node;
      % store basis function and gradients to speed up the code
      basis{p} = shape.phi;
      grad{p}  = shape.dphi;
      for i=1:length(esctr)
        id              = esctr(i);
        nmass(id)       = nmass(id)        + shape.phi(i)*mp;
        nmomentum0(id,:)= nmomentum0(id,:) + shape.phi(i)*mp*velo(p,:);
        niforce(id,1)   = niforce(id,1) - volp*(sig(1)*shape.dphi(i,1) + sig(3)*shape.dphi(i,2));
        niforce(id,2)   = niforce(id,2) - volp*(sig(3)*shape.dphi(i,1) + sig(2)*shape.dphi(i,2));
        neforce(id,:)   = neforce(id,:) + mp*shape.phi(i)*[bx, by]; 
      end
    end
    
    % debug
    % update nodal momenta    
    nforce    = niforce + neforce;
    nmomentum = nmomentum0 + nforce*dtime;
           
    % zero Dirichlet boundary conditions 
    
    nmomentum0(idx,:) = 0;
    nmomentum (idx,:) = 0;

    % update particle velocity and map back to nodes
    % loop over deformable particles
    for p=1:pCount   
        %data  = cpdi22D(p,particles,mesh); % old implementation
        input.nodes=nodeid{p};
        input.wf   = funcW{p};
        input.wg   = [gradW{p,1} gradW{p,2}];
        input.Vp   = volume(p);
  
        esctr   = input.nodes;
        phi     = basis{p};
        dphi    = grad{p};
        for i=1:length(esctr)
          id = esctr(i);
          if nmass(id) > 0
            velo(p,:)  = velo(p,:) + phi(i)*(nmomentum(id,:) - nmomentum0(id,:)) /nmass(id);            
          end          
        end
        % map back to nodes
        for i=1:length(esctr)
          id = esctr(i);          
          nvelo(id,:)  = nvelo(id,:) + phi(i)*velo(p,:)*mass(p);
        end
    end
    
    % Dirichlet boundary conditions
    nvelo(idx ,:) = 0;

    
    % loop over deformable particles
    for p=1:pCount
        Lp    = zeros(2,2);        
        %data  = cpdi22D(p,particles,mesh); % old implementation
        input.nodes=nodeid{p};
        input.wf   = funcW{p};
        input.wg   = [gradW{p,1} gradW{p,2}];
        input.Vp   = volume(p);
  
        esctr   = input.nodes;
        phi     = basis{p};
        dphi    = grad{p};
        for i=1:length(esctr)
          id = esctr(i);
          vI = [0 0];
          if nmass(id) > 0            
            vI = nvelo(id,:)/nmass(id);  % nodal velocity
          end
          Lp = Lp + vI'*dphi(i,:);         % particle gradient velocity
        end
        
        F          = ([1 0;0 1] + Lp*dtime)*reshape(deform(p,:),2,2);
        deform(p,:)= reshape(F,1,4);
        volume(p)  = det(F)*volume0(p);
        J          = det(F);
        b          = F*F';
        sigma        = 1/J*( mu*(b-I) + lambda*log(J)*I );
        stress(p,:)  = [sigma(1,1) sigma(2,2)  sigma(1,2) ];      
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
        if nmass(id) > 0
          xc = xc + dtime*N*nmomentum(id,:)/nmass(id);
        end
      end
      particles.node(c,:) = xc;
    end
    

    
    if (  mod(istep,intervalRemesh) == 0 )
      disp('negative Jacobian!!!Regenerating Voronoi cells')      
      [Node,Element,Supp,Load,P]=PolyMesher(@RingDomain,NElem,1,position);
      % store the particle mesh into a structure for convenience
      particles.node     = Node;
      particles.elem     = Element;
    end
    
        % update centroids (used to compute particle displacements)
    for e = 1:NElem
      vx  = particles.node(Element{e},1);
      vy  = particles.node(Element{e},2);
      nv  = length(Element{e});
      vxS = vx([2:nv 1]);
      vyS = vy([2:nv 1]); %Shifted vertices
      temp        = vx.*vyS - vy.*vxS;
      a           = 0.5*sum(temp);
      volume(e)   = a;
      position(e,:) = 1/(6*a)*[sum((vx+vxS).*temp),sum((vy+vyS).*temp)];
      u(e,:)        = position(e,:) - position0(e,:);
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
 
    % VTK output
    
    if (  mod(istep,interval) == 0 )
        vtkFile = sprintf('../../results/cpdi/%s%d',vtkFileName,istep);        
        data.stress  = [stress zeros(pCount,1)];        
        %data.disp    = u;
        data.color   = color;
        VTKPolygonParticles(particles,vtkFile,data);
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

% %
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
hold on
plot_mesh(node,element,'Q4','k-',1.); % background grid
patch('Faces',ElemMat,'Vertices',particles.node,'FaceColor','w'); pause(1e-6)
patch('Faces',ElemMat,'Vertices',Node,'FaceColor','b'); pause(1e-6)

disp([num2str(toc),'   DONE '])
