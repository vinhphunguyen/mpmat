% One dimensional bar in tension problem.
% Middle element is allowed  for cracking.
% Constitutive model: two-scale cohesive model with exponential softening.
% Dissipation based arc-length control to trace the snapback.
%
% VP Nguyen, vpnguyen@gmail.com, University of Adelaide, Australia.
% August, 2014.

addpath ../fem_util/
addpath ../gmshFiles/
addpath ../post-processing/
addpath ../fem-functions/

clear all
clc
state = 0;

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


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
mat.ft = 4;
mat.Gf = 20;
mat.ao = elasticityMatrix(E0,nu0,stressState);

option.implicit = 0;
option.tolerance= 1e-3;
option.tangent  = 1; %
option.iterMax  = 10;

% GENERATE FINITE ELEMENT MESH

lx       = 100;
ly       = 1;
numx     = 3;       % number of elements along X direction
numy     = 1;        % number of elements along Y direction
[mesh]   = buildGrid2D(lx,ly,numx,numy, 0);

element   = mesh.element;
node      = mesh.node;
nodeCount = mesh.nodeCount;
elemCount = mesh.elemCount;
mesh.sdof = 2*nodeCount;
mesh.elemType='Q4';
mesh.H    = mesh.deltax;%sqrt(mesh.deltax*mesh.deltay);

switchEnergy = 1e-5; % switch to arc-length indicator
switchIter   = 4;
optIter      = 3;
dEnergyMin   = 5e-4;
dEnergyMax   = 8e-3;


leftNodes  = find(node(:,1)==0);
botNodes   = find(node(:,2)==0);
rightNodes = find(abs(node(:,1)-lx)<1e-10);

% Essential boundary nodes (index and values)

udofs      = leftNodes;
vdofs      = botNodes+nodeCount;

uFixed     = zeros(length(leftNodes),1)';
vFixed     = zeros(length(botNodes),1)';

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

fhat(forcedDofs) = 1;

lambda0 = 0; % old force multiplier
dlambda = 0.04;
lambda  = dlambda; % current force multiplier

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
strain = zeros(3,4,elemCount); % strain vector (updated)
damage = zeros(  4,elemCount); % damage at Gauss points for visualisation

% ******************************************************************************
% ***                          P R O C E  S S I N G                          ***
% ******************************************************************************

noSteps  = 600;
tol      = 1e-5;
iterMax  = 60;  % maximum number of Newton-Raphson iterations

method   = 'forced-controlled';
diverged = 0;

rforce=[];
rdisp =[];

eStrain = [];
eStress = [];

used=0;

for i =1:3000
  
  disp ('=================================')
  disp (sprintf('   Load step %d',i));
  disp ('=================================')
  disp ('  NR iter : L2-norm residual')
  
  if     strcmp(method,'forced-controlled')
    disp (sprintf(' %s %%s %d ', ...
      'forced control with lambda: ', lambda) );
  else
    disp (sprintf(' %s %d %d ', ...
      'dissipation arc-length control with energy & lambda: ', dEnergy, lambda) );
  end
  diverged = 0;
  
  % compute the tangent and internal force
  
  [K,fint,status,ai,normal,jump,stress,strain] = ...
    computeTangentMatrix(mesh,mat,quad,u0,null,status,ai0,normal,jump0,stress0,option);
  
  error    = 1;
  iiter    = 0;
  
  du=zeros(mesh.sdof,1);
  
  while error > tol
    iiter    = iiter + 1;
    
    if     strcmp(method,'forced-controlled')
      % solve for the displacement increment
      res      = lambda*fhat-fint;
      [K,res]  = applyDirichletBCs(K,res,udofs,vdofs,uFixed,vFixed);
      ddu      = K\res;
    elseif strcmp(method,'energy-controlled')
      
      v     = 0.5*lambda0*fhat;
      omega = -0.5*u0'*fhat;
      phi   = 0.5*(lambda0*(u-u0)'-(lambda-lambda0)*u0')*fhat-dEnergy;
      %lambda,lambda0
      
      res      = lambda*fhat-fint;
      [K,res]  = applyDirichletBCs(K,res,udofs,vdofs,uFixed,vFixed);
      
      u1       = K\res;
      [K,fhat]  = applyDirichletBCs(K,fhat,udofs,vdofs,uFixed,vFixed);
      u2       = K\fhat;
      
      
      b        = v'*u1 + phi;
      c        = v'*u2 + omega;
      %
      %       if (c<1e-16)
      %         disp (sprintf(' %s %f ', 'c:', c) );
      %         error('Singular load increment!!!')
      %       end
      %phi
      %c
      alpha    = b/c;
      
      ddu      = u1 - alpha*u2;
      lambda   = lambda - alpha;
      %phi
    end
    
    % update displacements
    du = du + ddu;
    u  = u  + ddu;
    
    [K,fint,status,ai,normal,jump,stress,strain] = ...
      computeTangentMatrix(mesh,mat,quad,u0,du,status,ai0,normal,jump0,stress0,option);
    
    normf  = norm(lambda*fhat);
    
    error = norm( lambda*fhat-fint ) / normf;
    
    error = norm(ddu);
    
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
    %norm(u-u0)
    stress0 = stress;
    ai0     = ai;
    %ai0
    jump0   = jump;
    % check failure
    [status,normal] = checkFailure(mesh,mat,quad,stress0,status,ai0,normal);
    %ai0
    %status
    
    % store quantities for post-processing
    eStrain = [eStrain;strain(1,1,floor(mesh.elemCount/2)+1) ...
      strain(1,2,floor(mesh.elemCount/2)+1) ...
      strain(1,3,floor(mesh.elemCount/2)+1) ...
      strain(1,4,floor(mesh.elemCount/2)+1)];
    eStress = [eStress;stress(1,1,floor(mesh.elemCount/2)+1) ...
      stress(1,2,floor(mesh.elemCount/2)+1) ...
      stress(1,3,floor(mesh.elemCount/2)+1) ...
      stress(1,4,floor(mesh.elemCount/2)+1)];
    
    rforce = [rforce fint(forcedDofs(1))];
    rdisp  = [rdisp u0(forcedDofs(1))];
  end
  %   lambda0
  %   lambda
  dissnrg = 0.5 * ( lambda0 * du' - (lambda-lambda0) * u0'  )*fhat;
  
  disp (sprintf(' %s %d ','amount of dissipated energy: ', dissnrg));
  
  
  if     strcmp(method,'forced-controlled')
    if ~isempty(find(status==1, 1))
      disp('crack initiated!!!')
      %dlambda = dlambda/10;
    end
    if (dissnrg>switchEnergy) || (iiter >= switchIter) || ~isempty(find(status==1, 1))
      method='energy-controlled';
      dEnergy=dissnrg;
      %if (dEnergy<0.0001) dEnergy=0.0001;end
      if (dEnergy<dEnergyMin) dEnergy = dEnergyMin; end
      if (dEnergy>dEnergyMax) dEnergy = dEnergyMax; end
      disp (sprintf(' %s %d %s %d ', ...
        'switching to arc-length with dT=', dEnergy, ...
        'from load multiplier: ', lambda) );
      lambda0=lambda;
    else
      lambda0 = lambda;
      lambda  = lambda + dlambda;
    end
  elseif strcmp(method,'energy-controlled')
    lambda0 = lambda;
    dEnergy = dEnergy * 0.5^(0.25*(iiter-optIter));
    if (dEnergy<dEnergyMin) dEnergy = dEnergyMin; end
    if (dEnergy>dEnergyMax) dEnergy = dEnergyMax; end
  end
  
  
end

%save('bar160-25elems.mat','rforce','rdisp','eStrain','eStress');



% mesh sensitivity check

h1=load('bar100-15elems-good-Gf20ft4.mat');
h2=load('bar100-31elems-good-Gf20ft4.mat');

xx=500;
figure
set(gca,'FontSize',14)
hold on
plot(eStrain(1:end-xx,2),eStress(1:end-xx,2),'blue--','LineWidth',1.1);
%plot(h1.eStrain(1:end-xx,1),h1.eStress(1:end-xx,1),'black.','LineWidth',1.1);
%plot(h2.eStrain(1:end-xx,1),h2.eStress(1:end-xx,1),'black.','LineWidth',1.1);
% plot(eStrain(1:end-xx,2),eStress(1:end-xx,2),'blue--','LineWidth',1.1);
% plot(eStrain(1:end-xx,3),eStress(1:end-xx,3),'cy-.','LineWidth',1.6);
% plot(eStrain(1:end-xx,4),eStress(1:end-xx,4),'red--','LineWidth',1.6);
xlabel('strain')
ylabel('stress')
legend('GP1','GP2','GP3','GP4');
%axis([0 3 0 3])


figure
set(gca,'FontSize',14)
hold on
plot(rdisp(1:end-xx),rforce(1:end-xx)*2,'b--','LineWidth',1.6);
plot(h1.rdisp,h1.rforce*2,'black-','LineWidth',1.6);
plot(h2.rdisp,h2.rforce*2,'red-','LineWidth',1.6);
xlabel('displacement [mm]')
ylabel('force [N]')
legend('7 elements','15 elements','31 elements')
%axis([0 3 0 3])


figure
set(gca,'FontSize',14)
hold on
plot(eStrain(1:end-xx,2),eStress(1:end-xx,2),'blue--','LineWidth',1.1);
%plot(h1.eStrain(1:end-xx,1),h1.eStress(1:end-xx,1),'black-','LineWidth',1.1);
%plot(h2.eStrain(1:end-xx,1),h2.eStress(1:end-xx,1),'red-','LineWidth',1.1);
xlabel('strain')
ylabel('stress')
legend('7 elements','15 elements','31 elements')
%axis([0 3 0 3])

% scaleFact=1;
% xs=1:nodeCount;                  % x portion of u and v vectors
% ys=(nodeCount+1):2*nodeCount;      % y portion of u and v vectors
% Ux = u(xs);
% Uy = u(ys);
% plot_mesh(node+scaleFact*[Ux Uy],element,mesh.elemType,'g.-',1);


                                           % of element loop







disp([num2str(toc),'  RUN FINISHED'])
