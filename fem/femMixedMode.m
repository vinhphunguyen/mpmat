% Mixed-mode test Nooru-Mohamed.
% Constitutive model: two-scale cohesive model with exponential softening.
% Displacement control.
%
% Fast assembly using the triple format.
% Works for both triangular and quad elements of any order.
%
% VP Nguyen, vpnguyen@gmail.com, University of Adelaide, Australia.
% September, 2014.

addpath ../fem_util/
addpath ../fem-functions/
addpath ../post-processing/

clear all
clc
state = 0;

global crackedElems notchedElems cracks damage



% ******************************************************************************
% ***                            I N P  U T                                  ***
% ******************************************************************************
tic;
disp('************************************************')
disp('***          S T A R T I N G    R  U N        ***')
disp('************************************************')
disp([num2str(toc),'  START'])

% MATERIAL PROPERTIES

mat.E  = 32000;   % Young��s modulus
mat.nu = 0.2;  % Poisson ratio
mat.stressState ='PLANE_STRESS'; % set  to either 'PLANE_STRAIN' or "PLANE_STRESS'
mat.ft = 3.0;
mat.Gf = 0.11;
mat.ao = elasticityMatrix(mat.E,mat.nu,mat.stressState);
mat.penalty = 1e8; % penalty stiffness for cohesive law in compression

% OPTIONS FOR THE CONSTITUTIVE MODEL

option.implicit = 1;
option.tolerance= 1e-4;
option.tangent  = 1; %
option.iterMax  = 20;

% GENERATE FINITE ELEMENT MESH


meshFile = 'mixed-mode.msh';
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
mesh.thickness=50;

for e = 1:mesh.elemCount
    coord = mesh.node(mesh.element(e,:),:);
    a     = det([coord,[1;1;1]])/2;
    mesh.H(e) = sqrt(a);
end

% Finding node groups for boundary conditions

ngr1 = find(meshG.LINES(:,3)==1);
ngr2 = find(meshG.LINES(:,3)==2);
ngr3 = find(meshG.LINES(:,3)==3);

xyNodes = unique(meshG.LINES(ngr1,1:2));  % nodes fixed in both X and Y 
yfNodes  = unique(meshG.LINES(ngr2,1:2)); % nodes prescribed in y dir.
xfNodes  = unique(meshG.LINES(ngr3,1:2)); % nodes prescribed in x dir.

% Essential boundary nodes (index and values)

udofsFix   = xyNodes;
udofsFor   = xfNodes;
vdofsFix   = xyNodes+mesh.nodeCount;
vdofsFor   = yfNodes+mesh.nodeCount;
udofs      = [xyNodes;xfNodes];
vdofs      = [xyNodes;yfNodes]+mesh.nodeCount;

%PLOT MESH

figure
hold on
plot_mesh(mesh.node,mesh.element,mesh.elemType,'b-',1.1);
plot(mesh.node(xyNodes,1),mesh.node(xyNodes,2),'r*','markersize',15,'MarkerFaceColor','r');
plot(mesh.node(xfNodes,1),mesh.node(xfNodes,2),'bs','markersize',15,'MarkerFaceColor','b');
plot(mesh.node(yfNodes,1),mesh.node(yfNodes,2),'bo','markersize',15,'MarkerFaceColor','b');
axis off

% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),'  INITIALIZING DATA STRUCTURES'])

u0   = zeros (mesh.sdof,1);          % nodal displacement vector at old step
u    = zeros (mesh.sdof,1);          % nodal displacement vector (updated)
fext = zeros (mesh.sdof,1);          % external  load vector
fhat = zeros (mesh.sdof,1);          % unit external  load vector
null = zeros (mesh.sdof,1);

% for history variables

%[W,Q]=quadrature( 2, 'GAUSS', 2 ); % for quadrilaterals
[W,Q]=quadrature( 1, 'TRIANGULAR', 1 ); % for triangles, 1 GP for T3 elems

gpCount = length(W);

quad.W = W*mesh.thickness;
quad.Q = Q;

status   = zeros(  gpCount,mesh.elemCount); % cracked or not at all Gauss points
loading  = zeros(  gpCount,mesh.elemCount); % loading status at all Gauss points
loading0 = zeros(  gpCount,mesh.elemCount); % loading status at all Gauss points (previous converged step)
ai       = zeros(4,gpCount,mesh.elemCount); % tangent stiffness of localisation band
ai0      = zeros(4,gpCount,mesh.elemCount); % tangent stiffness of localisation band
normal   = zeros(2,gpCount,mesh.elemCount); % normal vectors at all Gauss points
kappa    = zeros(  gpCount,mesh.elemCount); % hostory variables (updated)
kappa0   = zeros(  gpCount,mesh.elemCount); % history variables previous converged step
jump     = zeros(2,gpCount,mesh.elemCount); % jump 
jump0    = zeros(2,gpCount,mesh.elemCount); % jump previous converged step
stress0  = zeros(3,gpCount,mesh.elemCount); % stress vector (old step)
stress   = zeros(3,gpCount,mesh.elemCount); % stress vector (updated)
strain   = zeros(3,gpCount,mesh.elemCount); % stress vector (updated)
traction = zeros(2,gpCount,mesh.elemCount); % traction vector (updated), for visualisation
damage   = zeros(  gpCount,mesh.elemCount); % damage variable for visualisation
cracks   = zeros(4,gpCount,mesh.elemCount);

% ******************************************************************************
% ***                          P R O C E  S S I N G                          ***
% ******************************************************************************


tol      = 1e-5;
iterMax  = 40;  % maximum number of Newton-Raphson iterations
da       = 0.0006;
ufinal   = 0.1;
noSteps  = floor(ufinal/da);

diverged = 0;

rforce=[];
rdisp =[];

eStrain = [];
eStress = [];

for i = 1:80
  
  %if (i>45) option.implicit = 1;end
  
  disp ('=================================')
  disp (sprintf('   Load step %d',i));
  disp ('=================================')
  disp (' NR iter: L2-norm residual')
  
  diverged = 0;
  
  % compute the tangent and internal force
  
  [K,fint,status,ai,normal,kappa,jump,stress,strain,traction,loading] = ...
    computeTangentMatrix(mesh,mat,quad,u0,null,status,ai0,normal,kappa0,jump0,stress0,loading0,option);
  
  error    = 1;
  iiter    = 0;
  
  du=zeros(mesh.sdof,1);
  
  while error > tol
    iiter    = iiter + 1;
    
    if (iiter==1)
      uFixed     = [zeros(length(udofsFix),1);-da*ones(length(udofsFor),1)]';
      vFixed     = [zeros(length(vdofsFix),1); da*ones(length(vdofsFor),1)]';
    else
      uFixed     = zeros(length(udofs),1)';
      vFixed     = zeros(length(vdofs),1)';
    end
    
    % solve for the displacement increment
    res      = fext-fint;
    [K,res]  = applyDirichletBCs(K,res,udofs,vdofs,uFixed,vFixed);
    ddu      = K\res;
    
    % update displacements
    du = du + ddu;
    u  = u + du;
    
    [K,fint,status,ai,normal,kappa,jump,stress,strain,traction,loading] = ...
      computeTangentMatrix(mesh,mat,quad,u0,du,status,ai0,normal,kappa0,jump0,stress0,loading0,option);
    
    if (iiter==1) 
      ddu1=ddu; 
      fint1=fint;
      rnmax=0;
    end
    
       
%     re    = fext-fint;
%     error = norm( re(freeDofs) );
%     F     = max ( norm(fext), norm(fint1) );
%     error = error/F;
%     
%     error = fint'*ddu/(fint1'*ddu1);
    
    rnorm = norm(ddu);
       rnmax = max(rnorm,rnmax);
    error = rnorm/norm(ddu1);
    
    disp (sprintf(' %s %i %s %5.4e ', 'Iter',iiter, ':', error) );
    
    if iiter == iterMax
      error('Newton-Raphson iterations did not converge!');
      diverged=1;
      break;
    end
  end
  
  %% the step has been converged...
  
  if (diverged==0)
    u0      = u0+du;
    stress0 = stress;
    kappa0  = kappa;
    jump0   = jump;
    ai0     = ai;
    loading0 = loading;
    % check failure
    [status,normal] = checkFailure(mesh,mat,quad,stress0,status,ai0,normal);
    
    for e=1:mesh.elemCount
      for q=1:size(quad.W,1)                          % quadrature loop
        if status (q,e)==0 continue;end
        %traction(1,q,e),kappa0(q,e),kappa(q,e)
        if (traction(1,q,e) < 0)
          %disp('ee')
          %status(q,e)=0;
        end
      end
    end
  end
  
  % store for load-displacement curve
  rforce = [rforce; sum(fint(vdofsFor)) sum(fint(udofsFor))];
  rdisp  = [rdisp; u0(vdofsFor(1)) u0(udofsFor(1))];
end

%save('mixed-.mat','','','cracks');

h1=load('beam867.mat');
h2=load('beam1663.mat');
h3=load('beam1663-no.mat');

h1.rdisp=[0 h1.rdisp];
h1.rforce=[0 h1.rforce];

load('mixed-mode-experiment.mat');
load('mixed-mode-experiment-shear.mat');
load('mixed-mode-giang.mat');
load('giang-shear.mat');

%plot load-displacement curve
xx=0;
figure
set(gca,'FontSize',14)
hold on
plot(rdisp(1:end-xx,1),0.001*rforce(1:end-xx,1),'red*-','LineWidth',1.6);
plot(exp(:,1),exp(:,2),'bo-','LineWidth',1.6);
plot(giang(:,1),giang(:,2),'r-','LineWidth',1.6);
xlabel('\delta_n [mm]')
ylabel('P_n [kN]')
legend('me','experiment','Giang')
%axis([0 0.5 0 1.2])

xx=0;
figure
set(gca,'FontSize',14)
hold on
plot(-rdisp(1:end-xx,2),-0.001*rforce(1:end-xx,2),'red*--','LineWidth',1.6);
plot(shear(:,1),shear(:,2),'bo-','LineWidth',1.6);
plot(unnamed(:,1),unnamed(:,2),'blacks-','LineWidth',1.6);
xlabel('\delta_s [mm]')
ylabel('P_s [kN]')
legend('me','experiment','Giang')
%axis([0 0.5 0 1.2])

%plot deformed mesh
figure
scaleFact=50;
xs=1:mesh.nodeCount;                  % x portion of u and v vectors
ys=(mesh.nodeCount+1):mesh.sdof;      % y portion of u and v vectors
Ux = u(xs);
Uy = u(ys);
plot_mesh(mesh.node+scaleFact*[Ux Uy],mesh.element,mesh.elemType,'black-',1);
hold on
plot_mesh(mesh.node+0*[Ux Uy],mesh.element,mesh.elemType,'b-',1);


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


% damage/stress at nodes by extrapolation
nodalDam   = zeros(mesh.elemCount,mesh.nodePerElem);
nodalSigXX = zeros(mesh.elemCount,mesh.nodePerElem);
nodalSigYY = zeros(mesh.elemCount,mesh.nodePerElem);
nodalSigXY = zeros(mesh.elemCount,mesh.nodePerElem);

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
sigmaVM = zeros(mesh.nodeCount,6);

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
VTKPostProcess(mesh.node,mesh.element,2,'Tri3','mixed-mode',...
    [sigmaVM(:,1) sigmaVM(:,2) sigmaVM(:,3) sigmaVM(:,4)],[Ux Uy]);

disp([num2str(toc),'  RUN FINISHED'])
