% Notched three point bending test.
% Constitutive model: two-scale cohesive model with exponential softening.
% Displacement control.
%
% Fast assembly using the triple format.
% Works for both triangular and quad elements of any order.
%
% VP Nguyen, vpnguyen@gmail.com, University of Adelaide, Australia.
% September, 2014.
%
% Status:
% 10/09/2014: displacement control with adaptive size.
%           : real time plot of load-displacement curve.
%           : write to VTK files.

addpath ../fem_util/
addpath ../fem-functions/
addpath ../post-processing/

clear all
clc
state = 0;

global crackedElems notchedElems cracks mat  a a0
global mesh quad u0 stress0 stress strain damage
global kappa  loading jump normal status
global kappa0 loading0 jump0

load('beam_experiment.mat'); % experiment load-displacement data
load('beam_experiment2.mat'); % experiment load-displacement data

vtkFileName='notched-beam-damage-local';

tic;

% ******************************************************************************
% ***                            I N P  U T                                  ***
% ******************************************************************************
tic;
disp('************************************************')
disp('***          S T A R T I N G    R  U N        ***')
disp('************************************************')
disp([num2str(toc),'  START'])

% MATERIAL PROPERTIES

mat.E  = 30000;   % Young��s modulus
mat.nu = 0.2;  % Poisson ratio
mat.stressState ='PLANE_STRESS'; % set  to either 'PLANE_STRAIN' or "PLANE_STRESS'
mat.ft0= 3.33;
mat.Gf = 0.124;
mat.ao = elasticityMatrix(mat.E,mat.nu,mat.stressState);
mat.penalty = 1e7; % penalty stiffness for cohesive law in compression
mat.ks      = 10;
mat.alpha = 0.99;
mat.beta  = 0.1;     % increading beta make the material brittle
mat.h     = mat.ft0^2/2/mat.Gf;

% OPTIONS FOR THE CONSTITUTIVE MODEL

option.implicit = 1;
option.tolerance= 1e-6;         % tolerance of implicit stress update
option.tangent  = 1;            %
option.iterMax  = 5;           % maximum # of NR iterations of implcit.
option.stepCount = 1;          % number of sub steps used in explicit update

% GENERATE FINITE ELEMENT MESH

meshFile = 'beam.msh';
meshG    = load_gmsh (meshFile);

elemType = 'T3';
mesh.nodeCount  = meshG.nbNod;
mesh.elemCount  = meshG.nbTriangles;
mesh.node     = meshG.POS(:,1:2);
mesh.element  = meshG.TRIANGLES(1:mesh.elemCount,1:3);
mesh.sdof     = 2*mesh.nodeCount;
mesh.elemType ='T3';
mesh.nodePerElem=3;
mesh.H        = zeros(mesh.elemCount,1);
mesh.thickness = 50;

for e = 1:mesh.elemCount
  coord = mesh.node(mesh.element(e,:),:);
  a     = det([coord,[1;1;1]])/2;
  mesh.H(e) = sqrt(a);
end

% Finding node groups for boundary conditions

ngr1 = find(meshG.POINTS(:,2)==1000);
ngr2 = find(meshG.POINTS(:,2)==2000);
ngr3 = find(meshG.POINTS(:,2)==3000);


fnodes = unique(meshG.POINTS(ngr1,1)); % nodes
node1  = unique(meshG.POINTS(ngr2,1)); % nodes
node2  = unique(meshG.POINTS(ngr3,1)); % nodes

% Essential boundary nodes (index and values)


udofs      = [node1];
vdofsFix   = [node1;node2]+mesh.nodeCount;
vdofsFor   = fnodes+mesh.nodeCount;
vdofs      = [node1;node2;fnodes]+mesh.nodeCount;

uFixed     = zeros(length(udofs),1)';

%PLOT MESH
plotMesh=0;
if plotMesh
  figure
  hold on
  plot_mesh(mesh.node,mesh.element,mesh.elemType,'black-',1.005);
  %plot(mesh.node(xyNodes,1),mesh.node(xyNodes,2),'r*','markersize',15,'MarkerFaceColor','r');
  %plot(mesh.node(xfNodes,1),mesh.node(xfNodes,2),'bs','markersize',15,'MarkerFaceColor','b');
  %plot(mesh.node(yfNodes,1),mesh.node(yfNodes,2),'bo','markersize',15,'MarkerFaceColor','b');
  axis off
end

% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),'  INITIALIZING DATA STRUCTURES'])

nsd  = 2;                            % number of spatial dimension
nstr = 3;                            % number of strain/stress component

u0   = zeros (mesh.sdof,1);          % nodal displacement vector at old step
u    = zeros (mesh.sdof,1);          % nodal displacement vector (updated)
fext = zeros (mesh.sdof,1);          % external  load vector
fhat = zeros (mesh.sdof,1);          % unit external  load vector
null = zeros (mesh.sdof,1);

xs=1:mesh.nodeCount;                  % x portion of u and v vectors
ys=(mesh.nodeCount+1):mesh.sdof;      % y portion of u and v vectors

% for history variables

%[W,Q]=quadrature( 2, 'GAUSS', 2 ); % for quadrilaterals
[W,Q]=quadrature( 1, 'TRIANGULAR', 1 ); % for triangles, 1 GP for T3 elems

gpCount = length(W);

quad.W = W*mesh.thickness;
quad.Q = Q;

status   = zeros(       gpCount,mesh.elemCount); % cracked or not at all Gauss points
loading  = zeros(       gpCount,mesh.elemCount); % loading status at all Gauss points
loading0 = zeros(       gpCount,mesh.elemCount); % loading status at all Gauss points (previous converged step)
ai       = zeros(9, gpCount,mesh.elemCount); % tangent stiffness of localisation band
ai0      = zeros(9, gpCount,mesh.elemCount); % tangent stiffness of localisation band
a        = zeros(9,gpCount,mesh.elemCount); % homogenised tangent stiffness
a0       = zeros(9,gpCount,mesh.elemCount); % homogenised tangent stiffness
normal   = zeros(nsd,   gpCount,mesh.elemCount); % normal vectors at all Gauss points
kappa    = zeros(       gpCount,mesh.elemCount); % history variables (updated)
kappa0   = zeros(       gpCount,mesh.elemCount); % history variables previous converged step
jump     = zeros(3,   gpCount,mesh.elemCount); % jump
jump0    = zeros(3,   gpCount,mesh.elemCount); % jump previous converged step
stress0  = zeros(nstr,  gpCount,mesh.elemCount); % stress vector (old step)
stress   = zeros(nstr,  gpCount,mesh.elemCount); % stress vector (updated)
strain   = zeros(nstr,  gpCount,mesh.elemCount); % stress vector (updated)
traction = zeros(nsd,   gpCount,mesh.elemCount); % traction vector (updated), for visualisation
damage   = zeros(       gpCount,mesh.elemCount); % damage variable for visualisation
mat.fts  = zeros(       gpCount,mesh.elemCount); % variable tensile strength to
% ensure traction contuinuity
% at crack initiation
mat.kis  = (mat.ft0/mat.E)*ones (  gpCount,mesh.elemCount);

nodalDam   = zeros(mesh.elemCount,mesh.nodePerElem);
nodalSigXX = zeros(mesh.elemCount,mesh.nodePerElem);
nodalSigYY = zeros(mesh.elemCount,mesh.nodePerElem);
nodalSigXY = zeros(mesh.elemCount,mesh.nodePerElem);
sigmaVM    = zeros(mesh.nodeCount,6);

% ******************************************************************************
% ***                          P R O C E  S S I N G                          ***
% ******************************************************************************

disp([num2str(toc),'  SOLVING...'])

tol      = 1e-6;
iterMax  = 15;  % maximum number of Newton-Raphson iterations
da       = 0.01;
ufinal   = 0.8;
daMin    = 2e-4;
daMax    = 1e-2;

noSteps  = floor(ufinal/da);

resolve  = 1; % resolve the load increment with new cracks
diverged = 0;
optIter  = 7;

rforce=[];
rdisp =[];

eStrain = [];
eStress = [];

% real time plot of load-displacement curve
figure
clf;
set( gcf, 'DoubleBuffer', 'on' )
set(gca,'FontSize',14)
hold on
plot(beam_experiment(:,1),beam_experiment(:,2),'bo-','LineWidth',1.6);
plot(beam_exp_2(:,1),beam_exp_2(:,2),'bo-','LineWidth',1.6);
xlabel('u [mm]')
ylabel('P [N]')
%legend('fixed, 867 elems','fixed, 1633 elems','no fixed, 1633 elems')
xlim([0,0.8]); ylim([0,800]);
drawnow
%linkdata on

ii = 1;
isResolved=0;
interval = 4;
ipara=0;
for i = 1:noSteps*100
  
  disp ('=================================')
  disp (sprintf('   Load step %d',ii));
  disp ('=================================')
  disp (' NR iter: L2-norm residual')
  
  diverged = 0;
  
  % compute the tangent and internal force
  [K,fint,ai,traction] = computeTangentMatrix(null,ai0,option);
  
  
  error    = 1;
  iiter    = 0;
  
  du=zeros(mesh.sdof,1);
  
  while (error > tol) && (iiter<iterMax)
    iiter    = iiter + 1;
    
    if (iiter==1)
      vFixed     = [zeros(length(vdofsFix),1); -da*ones(length(vdofsFor),1)]';
    else
      vFixed     = zeros(length(vdofs),1)';
    end
    
    % solve for the displacement increment
    res      = fext-fint;
    [K,res]  = applyDirichletBCs(K,res,udofs,vdofs,uFixed,vFixed);
    ddu      = K\res;
    
    % update displacements
    du = du + ddu;
    u  = u + du;
    
    [K,fint,ai,traction] = computeTangentMatrix(du,ai0,option);
    
    if (iiter==1)
      ddu1=ddu;
      fint1=fint;
      rnmax=0;
    end
    
    rnorm = norm(ddu);
    rnmax = max(rnorm,rnmax);
    error = rnorm/rnmax;
    
    disp (sprintf(' %s %i %s %5.4e ', 'Iter',iiter, ':', error) );
    
    if iiter == iterMax
      warning('Newton-Raphson iterations did not converge!');
      diverged=1;
      %break;
    end
  end
  
  %% the step has been converged...
  if (diverged==0)
    % check failure
    
    isFailure=0;
    %if (resolve==0) isFailure=0; end
    
    % if no cracks is initiated then proceed to next load increment
    if (isFailure==0)
      u0       = u0+du;
      stress0  = stress;
      kappa0   = kappa;
      jump0    = jump;
      loading0 = loading;
      ai0      = ai;
      a0       = a;
      ii       = ii + 1;
      
      % store for load-displacement curve
      rforce = [rforce; fint(vdofsFor(1)) ];
      rdisp  = [rdisp; u0(vdofsFor(1)) ];
      
      % real time plot of the load-displacement curve
      plot(-rdisp(1:end),-rforce(1:end),'black-*','LineWidth',1.6);
      xlim([0,0.8]); ylim([0,800]);
      drawnow
      
      if (  mod(ii-1,interval) == 0 )
        % damage/stress at nodes by extrapolation
        nodalDam(:)   = 0;
        nodalSigXX(:) = 0;
        nodalSigYY(:) = 0;
        nodalSigXY(:) = 0;
        
        for e=1:mesh.elemCount
          for q=1:size(quad.W,1)                          % quadrature loop
            pt       = Q(q,:);                       % quadrature point
            wt       = W(q);                         % quadrature weight
            [N,~]         = lagrange_basis(mesh.elemType,pt);
            d             = damage(q,e);
            nodalDam(e,:)   =  nodalDam(e,:)   + N'*d;
            nodalSigXX(e,:) =  nodalSigXX(e,:) + N'*stress(1,q,e);
            nodalSigYY(e,:) =  nodalSigYY(e,:) + N'*stress(2,q,e);
            nodalSigXY(e,:) =  nodalSigXY(e,:) + N'*stress(3,q,e);
          end                                               % of element loop
        end
        
        % avaraging nodal damage
        sigmaVM(:) = 0;
        for e=1:mesh.elemCount
          connect = mesh.element(e,:);
          for in=1:mesh.nodePerElem
            nid = connect(in);
            sigmaVM(nid,1)     = sigmaVM(nid,1) + nodalDam  (e,in);
            sigmaVM(nid,2)     = sigmaVM(nid,2) + nodalSigXX(e,in);
            sigmaVM(nid,3)     = sigmaVM(nid,3) + nodalSigYY(e,in);
            sigmaVM(nid,4)     = sigmaVM(nid,4) + nodalSigXY(e,in);
            sigmaVM(nid,5)     = sigmaVM(nid,5) + 1;
            if nodalDam(e,in)~=0
              sigmaVM(nid,6) = sigmaVM(nid,6) + 1;
            end
          end
        end
        
        % Average nodal stress values (learned from Mathiew Pais XFEM code)
        
        for i=1:mesh.nodeCount
          if sigmaVM(i,6)~=0
            sigmaVM(i,1) = sigmaVM(i,1)/sigmaVM(i,6);
          end
        end
        
        sigmaVM(:,2) = sigmaVM(:,2)./sigmaVM(:,5);
        sigmaVM(:,3) = sigmaVM(:,3)./sigmaVM(:,5);
        sigmaVM(:,4) = sigmaVM(:,4)./sigmaVM(:,5);
        
        %sigmaVM(:,5) = [];
        %sigmaVM(:,6) = [];
        
        % write to Paraview
        Ux = u0(xs);
        Uy = u0(ys);
        
        vtkFile = sprintf('%s%d',vtkFileName,ipara);
        ipara=ipara+1;
        VTKPostProcess(mesh.node,mesh.element,2,'Tri3',vtkFile,...
          [sigmaVM(:,1) sigmaVM(:,2) sigmaVM(:,3) sigmaVM(:,4)],[Ux Uy]);
      end                                 % end of writing VTK
      %da = da * 0.5^(0.25*(iiter-optIter));
    end                                   % end of isFailure==0 check
  else                                    % adaptive reducing step size
    da = da * 0.5^(0.25*(iiter-optIter));
    
    if (da < daMin), da=daMin; end
    if (da > daMax), da=daMax; end
  end                                     % end of divergence check
  
end

disp([num2str(toc),'  SOLVING FINISHED'])

%%

disp([num2str(toc),'  POST-PROCESSING'])

%%


pvdFile = fopen([vtkFileName '.pvd'], 'wt');

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');

for i = 1:ipara
  vtuFile = sprintf('%s%d%s',vtkFileName,i,'.vtu');
  fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
end

fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);

%%

Ux = u0(xs);
Uy = u0(ys);

%plot load-displacement curve
xx=0;
figure
set(gca,'FontSize',14)
hold on
plot(-rdisp(1:end-1),-rforce(1:end-1),'red*-','LineWidth',1.6);
plot(beam_experiment(:,1),beam_experiment(:,2),'bo-','LineWidth',1.6);
plot(beam_exp_2(:,1),beam_exp_2(:,2),'bo-','LineWidth',1.6);
xlabel('Deflection [mm]')
ylabel('Load P [N]')
legend('me','experiment')
axis([0 0.8 0 800])


%plot deformed mesh
figure
scaleFact=50;
plot_mesh(mesh.node+scaleFact*[Ux Uy],mesh.element,mesh.elemType,'black-',1);
hold on
%plot_mesh(mesh.node+0*[Ux Uy],mesh.element,mesh.elemType,'b-',1);


%plot undeformed mesh with active cracks (to see crack pattern)
figure
hold on
plot_mesh(mesh.node,mesh.element,mesh.elemType,'b-',1.);
for e=1:mesh.elemCount
  for q=1:size(quad.W,1)                          % quadrature loop
    if status (q,e)==0 continue; end
    if loading(q,e)==0 continue; end
    plot(cracks([1 3],q,e),cracks([2 4],q,e),'black-','LineWidth',1.8)
  end
end
axis off

%plot undeformed mesh with both active and closed cracks (to see crack pattern)
figure
hold on
plot_mesh(mesh.node,mesh.element,mesh.elemType,'b-',1.);
for e=1:mesh.elemCount
  for q=1:size(quad.W,1)                          % quadrature loop
    if status (q,e)==0 continue; end
    plot(cracks([1 3],q,e),cracks([2 4],q,e),'red--','LineWidth',1.8)
  end
end
axis off




disp([num2str(toc),'  RUN FINISHED'])
