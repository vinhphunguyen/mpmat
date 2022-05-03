% Notched three point bending problem with Q4 elements.
% Constitutive model: two-scale cohesive model with exponential softening.
% Displacement control to trace the snapback.
%
% Fast assembly using the triple format.
% Works for both triangular and quad elements of any order.
%
% VP Nguyen, vpnguyen@gmail.com, University of Adelaide, Australia.
% 13/September, 2014.

addpath ../fem_util/
addpath ../fem-functions/
addpath ../post-processing/

clear all
clc
state = 0;


opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

global crackedElems notchedElems cracks mat  a a0
global mesh quad u0 stress0 stress strain damage
global kappa  loading jump normal status
global kappa0 loading0 jump0

cracks=[];

tic;

load('beam_experiment.mat'); % experiment load-displacement data
load('beam_experiment2.mat'); % experiment load-displacement data

% ******************************************************************************
% ***                            I N P  U T                                  ***
% ******************************************************************************
tic;
disp('************************************************')
disp('***          S T A R T I N G    R  U N        ***')
disp('************************************************')
disp([num2str(toc),'  START'])

% MATERIAL PROPERTIES

mat.E  = 30000;   % Young€™s modulus
mat.nu = 0.2;  % Poisson ratio
mat.stressState ='PLANE_STRESS'; % set  to either 'PLANE_STRAIN' or "PLANE_STRESS'
mat.ft0= 3.33;
mat.Gf = 0.124;
mat.ao = elasticityMatrix(mat.E,mat.nu,mat.stressState);
mat.penalty = 1e7; % penalty stiffness for cohesive law in compression
mat.ks      = 10;

% OPTIONS FOR THE CONSTITUTIVE MODEL

option.implicit = 1;
option.tolerance= 1e-7;
option.tangent  = 1; %
option.iterMax  = 10;
option.stepCount = 1;

vtkFileName='unnotched-beam-explicit';

% GENERATE FINITE ELEMENT MESH

thickness=1;
lx       = 2000;
ly       = 200;
numx     = 201;       % number of elements along X direction
numy     = 20;       % number of elements along Y direction
[mesh]   = buildGrid2D(lx,ly,numx,numy, 0);

mesh.sdof = 2*mesh.nodeCount;
mesh.elemType='Q4';
mesh.H    = sqrt(mesh.deltax*mesh.deltay)*ones(mesh.elemCount,1);
mesh.nodePerElem=4;
mesh.thickness = 50;

node1      = 1;
node2      = numx + 1;

x1 = mesh.deltax*floor(numx/2);
x2 = x1 + mesh.deltax;
tt1=find(abs(mesh.node(:,1)-x1)<1e-10);
tt2=find(abs(mesh.node(:,1)-x2)<1e-10);

fnodes=[tt1(end);tt2(end)];

xNodes=[node1];
yNodes=[node1;node2;fnodes];
xDofs = xNodes;
yDofs = yNodes+mesh.nodeCount;
consDofs=[xDofs;yDofs];
freeDofs=setdiff(1:mesh.sdof,consDofs);

% force cracking in the column of elements
crackedElems=floor(numx/2)+1:numx:mesh.elemCount;
notchedElems=crackedElems(1:numy/2);
crackedElems=[];
% cracking everywhere except at the support
%crackedElems=setdiff(1:mesh.elemCount,[1 2 3 4 5 numx-3 numx-4 numx-2 numx-1 numx]);

% Essential boundary nodes (index and values)

udofs      = [node1];
vdofsFix   = [node1;node2]+mesh.nodeCount;
vdofsFor   = fnodes+mesh.nodeCount;
vdofs      = [node1;node2;fnodes]+mesh.nodeCount;

uFixed     = zeros(length(udofs),1)';

%PLOT MESH

figure
hold on
plot_mesh(mesh.node,mesh.element,mesh.elemType,'b-',1.);
plot_mesh(mesh.node,mesh.element(notchedElems,:),mesh.elemType,'cy-',1.1);
plot_mesh(mesh.node,mesh.element(crackedElems,:),mesh.elemType,'r-',1.1);
% plot(node(leftNodes1,1),node(leftNodes1,2),'s');
% plot(node(node2,1),node(node2,2),'s');
axis off

% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),'  INITIALIZING DATA STRUCTURES'])

u0   = zeros (mesh.sdof,1);          % nodal displacement vector at old step
u    = zeros (mesh.sdof,1);          % nodal displacement vector (updated)
fext = zeros (mesh.sdof,1);          % external  load vector
fhat = zeros (mesh.sdof,1);          % unit external  load vector
null = zeros (mesh.sdof,1);

xs=1:mesh.nodeCount;                  % x portion of u and v vectors
ys=(mesh.nodeCount+1):mesh.sdof;      % y portion of u and v vectors

% for history variables

[W,Q]=quadrature( 2, 'GAUSS', 2 ); % for quadrilaterals
%[W,Q]=quadrature( 1, 'TRIANGULAR', 1 ); % for triangles, 1 GP for T3 elems

gpCount = length(W);

quad.W = W*mesh.thickness;
quad.Q = Q;

status   = zeros(  gpCount,mesh.elemCount); % cracked or not at all Gauss points
loading  = zeros(  gpCount,mesh.elemCount); % loading status at all Gauss points
loading0 = zeros(  gpCount,mesh.elemCount); % loading status at all Gauss points (previous converged step)
ai       = zeros(4,gpCount,mesh.elemCount); % tangent stiffness of localisation band
ai0      = zeros(4,gpCount,mesh.elemCount); % tangent stiffness of localisation band
a        = zeros(9,gpCount,mesh.elemCount); % homogenised tangent stiffness
a0       = zeros(9,gpCount,mesh.elemCount); % homogenised tangent stiffness
normal   = zeros(2,gpCount,mesh.elemCount); % normal vectors at all Gauss points
kappa    = zeros(  gpCount,mesh.elemCount); % history variables (updated)
kappa0   = zeros(  gpCount,mesh.elemCount); % history variables previous converged step
jump     = zeros(2,gpCount,mesh.elemCount); % jump
jump0    = zeros(2,gpCount,mesh.elemCount); % jump previous converged step
stress0  = zeros(3,gpCount,mesh.elemCount); % stress vector (old step)
stress   = zeros(3,gpCount,mesh.elemCount); % stress vector (updated)
strain   = zeros(3,gpCount,mesh.elemCount); % stress vector (updated)
traction = zeros(2,gpCount,mesh.elemCount); % traction vector (updated), for visualisation
damage   = zeros(  gpCount,mesh.elemCount); % damage variable for visualisation
cracks   = zeros(4,gpCount,mesh.elemCount);
mat.fts  = zeros(  gpCount,mesh.elemCount); % variable tensile strength to
% ensure traction contuinuity
% at crack initiation

nodalDam   = zeros(mesh.elemCount,mesh.nodePerElem);
nodalSigXX = zeros(mesh.elemCount,mesh.nodePerElem);
nodalSigYY = zeros(mesh.elemCount,mesh.nodePerElem);
nodalSigXY = zeros(mesh.elemCount,mesh.nodePerElem);
sigmaVM    = zeros(mesh.nodeCount,6);

% ******************************************************************************
% ***                          P R O C E  S S I N G                          ***
% ******************************************************************************


tol      = 1e-8;
iterMax  = 15;  % maximum number of Newton-Raphson iterations
da       = 0.01;
ufinal   = 0.8;
noSteps  = floor(ufinal/da);
daMin    = 2e-4;
daMax    = 1e-2;
optIter  = 7;

interval = 2;

rforce=[0];
rdisp =[0];

eStrain = [];
eStress = [];

ii = 1;             % load increment index, there are steps that are resolved
isResolved=0;
ipara=0;

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

for i = 1:noSteps
  
  disp ('=================================')
  disp (sprintf('   Load step %d',ii));
  disp ('=================================')
  disp ('  NR iter : L2-norm residual')
  
  diverged = 0;
  
  % compute the tangent and internal force
  [K,fint,ai,traction] = computeTangentMatrix(null,ai0,option);
  
  error    = 1;
  iiter    = 0;
  
  du=zeros(mesh.sdof,1);
  
  while error > tol
    iiter    = iiter + 1;
    
    if (iiter==1)
      vFixed     = [zeros(length(vdofsFix),1);-da*ones(length(vdofsFor),1)]';
    else
      vFixed     = zeros(length(vdofs),1)';
    end
    
    % solve for the displacement increment
    res      = fext-fint;
    [K,res]  = applyDirichletBCs(K,res,udofs,vdofs,uFixed,vFixed);
    ddu      = K\res;
    
    % update displacements
    du = du + ddu;
    u  = u  + du;
    
    [K,fint,ai,traction] = computeTangentMatrix(du,ai0,option);
    
    if (iiter==1)
      ddu1=ddu;
      fint1=fint;
      rnmax=0;
    end
    
    % displacement norm
    rnorm = norm(ddu);
    rnmax = max(rnorm,rnmax);
    error = rnorm/rnmax;
    
    % force norm
    %error = norm(res)/norm(fint1);
    
    disp (sprintf(' %s %i %s %5.4e ', 'Iter',iiter, ':', error) );
    
    if iiter == iterMax
      warning('Newton-Raphson iterations did not converge!');
      diverged=1;%da=da/2;
      break;
    end
  end
  
  %% the step has been converged...
  if (diverged==0)
    % check failure
    isFailure = checkFailure(ai0);
    %isFailure=0;
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
      
      % switch compressed crack points to continuum elastic model
      %       for e=1:mesh.elemCount
      %         for q=1:size(quad.W,1)                          % quadrature loop
      %           if status (q,e)==0 continue;end
      %           %traction(1,q,e),kappa0(q,e),kappa(q,e)
      %           if (traction(1,q,e) < 0)
      %             %disp('ee')
      %             %status(q,e)=0;
      %           end
      %         end
      %       end
      
      % store for load-displacement curve
      rforce = [rforce; fint(vdofsFor(1)) ];
      rdisp  = [rdisp; u0(vdofsFor(1)) ];
      
      % real time plot of the load-displacement curve
      plot(-rdisp(1:end),-rforce(1:end)*2,'black-*','LineWidth',1.);
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
        Ux = u0(xs);
        Uy = u0(ys);
        % write to Paraview
        vtkFile = sprintf('%s%d',vtkFileName,ipara);
        ipara=ipara+1;
        VTKPostProcess(mesh.node,mesh.element,2,'Quad4',vtkFile,...
          [sigmaVM(:,1) sigmaVM(:,2) sigmaVM(:,3) sigmaVM(:,4)],[Ux Uy]);
      end
      %isResolved=0;
      %else
      %isResolved=1;
    end
  else
    da = da * 0.5^(0.25*(iiter-optIter));
  end
end


disp([num2str(toc),'  SOLVING FINISHED'])

%%

disp([num2str(toc),'  POST-PROCESSING'])

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
figure
clf;
set( gcf, 'DoubleBuffer', 'on' )
set(gca,'FontSize',14)
hold on
plot(beam_experiment(:,1),beam_experiment(:,2),'bo-','LineWidth',1.6);
plot(beam_exp_2(:,1),beam_exp_2(:,2),'bo-','LineWidth',1.6);
plot(-rdisp(1:end),-rforce(1:end)*2,'black-*','LineWidth',1.);
xlabel('u [mm]')
ylabel('P [N]')
%legend('fixed, 867 elems','fixed, 1633 elems','no fixed, 1633 elems')
xlim([0,0.8]); ylim([0,800]);
drawnow
%linkdata on

%plot deformed mesh
figure
scaleFact=1000;
plot_mesh(mesh.node+scaleFact*[Ux Uy],mesh.element,mesh.elemType,'black-',1);

%plot undeformed mesh with active cracks (to see crack pattern)
figure
hold on
plot_mesh(mesh.node,mesh.element,mesh.elemType,'b-',1.);
for e=1:mesh.elemCount
  for q=1:size(quad.W,1)                          % quadrature loop
    if status (q,e)==0 continue; end
    if loading0(q,e)==0 continue; end
    plot(cracks([1 3],q,e),cracks([2 4],q,e),'black-','LineWidth',1.8)
  end
end
axis off

%plot undeformed mesh with both active and closed cracks (to see crack pattern)
figure
hold on
plot_mesh(mesh.node,mesh.element,mesh.elemType,'b-',1.);
plot_mesh(mesh.node,mesh.element(notchedElems,:),mesh.elemType,'cy-',1.1);
for e=1:mesh.elemCount
  for q=1:size(quad.W,1)                          % quadrature loop
    if status (q,e)==0 continue; end
    plot(cracks([1 3],q,e),cracks([2 4],q,e),'red--','LineWidth',1.8)
  end
end
axis off

disp([num2str(toc),'  RUN FINISHED'])
