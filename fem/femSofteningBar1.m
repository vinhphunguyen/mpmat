

addpath ../fem_util/
addpath ../fem-functions/
addpath ../../post-processing/

clear all
clc
state = 0;

% ******************************************************************************
% ***                            I N P  U T                                  ***
% ******************************************************************************
tic;
disp('************************************************')
disp('***          S T A R T I N G    R  U N        ***')
disp('************************************************')
disp([num2str(toc),'  START'])

% MATERIAL PROPERTIES

E0  = 10;   % Young€™s modulus
nu0 = 0.0;  % Poisson ratio
stressState ='PLANE_STRAIN'; % set  to either 'PLANE_STRAIN' or "PLANE_STRESS'
mat.ft = 1;
mat.Gf = 1;
mat.ao = elasticityMatrix(E0,nu0,stressState);


option.implicit = 1;
option.tolerance= 0.001;
option.tangent  = 1; %

% GENERATE FINITE ELEMENT MESH

lx       = 4;
ly       = 1;
numx     = 17;       % number of elements along X direction
numy     = 1;        % number of elements along Y direction
[mesh]   = buildGrid2D(lx,ly,numx,numy, 0);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;
mesh.sdof = 2*nodeCount;
mesh.elemType='Q4';
mesh.H    = mesh.deltax;%sqrt(mesh.deltax*mesh.deltay);


leftNodes  = find(node(:,1)==0);
rightNodes = find(abs(node(:,1)-lx)<1e-10);

% Essential boundary nodes (index and values)

udofs      = [leftNodes;rightNodes];
vdofs      = leftNodes+nodeCount;

vFixed     = zeros(length(leftNodes),1)';

forcedDofs = rightNodes;

%PLOT MESH

clf
plot_mesh(node,element,mesh.elemType,'g.-',1.1);
hold on
axis off

% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),'  INITIALIZING DATA STRUCTURES'])

u0   = zeros (mesh.sdof,1);          % nodal displacement vector at old step
u    = zeros (mesh.sdof,1);          % nodal displacement vector (updated)
fext = zeros (mesh.sdof,1);          % external  load vector
fhat = zeros (mesh.sdof,1);          % unit external  load vector
null = zeros (mesh.sdof,1);

% for history variables

[W,Q]=quadrature( 2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature
quad.W = W;
quad.Q = Q;

status = zeros(  4,elemCount); % cracked or not at all Gauss points
ai     = zeros(4,4,elemCount); % tangent stiffness of localisation band
ai0    = zeros(4,4,elemCount); % tangent stiffness of localisation band
normal = zeros(2,4,elemCount); % normal vectors at all Gauss points
jump   = zeros(2,4,elemCount); % displacement jump at all Gauss points
jump0  = zeros(2,4,elemCount); % displacement jump at all Gauss points
stress0= zeros(3,4,elemCount); % stress vector (old step)
stress = zeros(3,4,elemCount); % stress vector (updated)
strain = zeros(3,4,elemCount); % stress vector (updated)

% ******************************************************************************
% ***                          P R O C E  S S I N G                          ***
% ******************************************************************************


tol      = 1e-6;
iterMax  = 30;  % maximum number of Newton-Raphson iterations
da       = 0.006;
ufinal   = 1*lx/2;
noSteps  = floor(ufinal/da);

diverged = 0;

rforce=[];
rdisp =[];

eStrain = [];
eStress = [];

for i = 1:noSteps
  
  disp ('=================================')
  disp (sprintf('   Load step %d',i));
  disp ('=================================')
  disp ('  NR iter : L2-norm residual')
  
  diverged = 0;
  
  % compute the tangent and internal force
  
  [K,fint,status,ai,normal,jump,stress,strain] = ...
    computeTangentMatrix(mesh,mat,quad,u0,null,status,ai,normal,jump0,stress0,option);
  
  error    = 1;
  iiter    = 0;
  
  du=zeros(mesh.sdof,1);
  
  while error > tol
    iiter    = iiter + 1;
    
    if (iiter==1)
      uFixed     = [zeros(length(leftNodes),1);da*ones(length(rightNodes),1)]';
    else
      uFixed     = zeros(length(udofs),1)';
    end
    
    % solve for the displacement increment
    res      = fext-fint;
    [K,res]  = applyDirichletBCs(K,res,udofs,vdofs,uFixed,vFixed);
    ddu      = K\res;
    
    % update displacements
    du = du + ddu;
    u  = u + du;
    
    [K,fint,status,ai,normal,jump,stress,strain] = ...
      computeTangentMatrix(mesh,mat,quad,u0,du,status,ai,normal,jump0,stress0,option);
    
    normf  = norm(fext);
    if normf < 1e-16
      error = norm( fext-fint );
    else
      error = norm( fext-fint ) / normf;
    end
    
    if (iiter==1) ddu1=ddu; end
    
    rnorm = norm(ddu);
%    rnmax = max(rnorm,rnmax);
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
    jump0   = jump;
    % check failure
    [status,ai,normal] = checkFailure(mesh,mat,quad,stress0,status,ai,normal);
    
    ai;
  end
  
  eStress = [eStress;stress(1,1,floor(mesh.elemCount/2)+1) stress(1,2,floor(mesh.elemCount/2)+1)];
  eStrain = [eStrain;strain(1,1,floor(mesh.elemCount/2)+1) strain(1,2,floor(mesh.elemCount/2)+1)];
  
  rforce = [rforce fint(forcedDofs(1))];
  rdisp  = [rdisp u0(forcedDofs(1))];
end

figure
set(gca,'FontSize',14)
hold on
plot(rdisp,2*rforce,'black-','LineWidth',1.6);
xlabel('strain')
ylabel('stress')
legend('\alpha=1','\alpha=100','\alpha=1000','\alpha=1000,explicit')
%axis([0 3 0 3])

figure
set(gca,'FontSize',14)
hold on
plot(eStrain(1:end-0,1),eStress(1:end-0,1),'black-','LineWidth',1.6);
plot(eStrain(1:end-0,2),eStress(1:end-0,2),'blue*-','LineWidth',1.6);
xlabel('strain')
ylabel('stress')
legend('GP1','GP2');
%axis([0 3 0 3])


%save('barElem3.mat','rdisp','rforce');
h3=load('barElem3.mat');
h7=load('barElem7.mat');
h15=load('barElem15.mat');

figure
set(gca,'FontSize',14)
hold on
plot(h3.rdisp,h3.rforce*2,'cyan*','LineWidth',1.6);
plot(h7.rdisp,h7.rforce*2,'black-','LineWidth',1.6);
plot(h15.rdisp,h15.rforce*2,'r-','LineWidth',1.6);
plot(rdisp,rforce*2,'b-','LineWidth',1.6);
xlabel('strain')
ylabel('stress')
legend('3 elemts','7 elems','15 elems','implcit')
%axis([0 3 0 3])

disp([num2str(toc),'  RUN FINISHED'])
