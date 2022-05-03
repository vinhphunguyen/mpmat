% An example of the two-scale cohesive model.
% Implicit stress update.
% Bulk: elastic
% Crack: initially rigid exponential cohesive law.
% VP Nguyen, phu.nguyen@adelaide.edu.au
% University of Adelaide, Australia, August, 2014.

addpath ../fem_util/
addpath ../fem-functions/
addpath ../post-processing/

clc;

% used to export fig to EPS file
opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

global mat

mat.E  = 10;
mat.nu = 0.;
mat.ft = 1;
mat.Gf = 1;
mat.ao = elasticityMatrix(mat.E,mat.nu,'PLANE_STRAIN');
mat.ks = 0;
mat.penalty=1e3;
mat.stressState ='PLANE_STRESS';

H  = 4;

% option for the stress update algorithm
option.implicit = 1;
option.tolerance= 1e-8;
option.tangent  = 1; %
option.iterMax  = 15;
option.stepCount = 1;

% initial values
status   = 0;
loading0 = 0;
kappa0   = 0;
jump0    = [0;0];
stress0  = [0;0;0];
strain0  = [0;0;0];
a0       = zeros(9,1);
ai       = zeros(4,1);
normal   = [1;0];

% loadings
epsilonF = mat.ft/mat.E;
epsilonf = 1;
epsilonr = 0.2;
dEpsE    = 1e-4;
alpha    = 100;
dEpsN    = dEpsE*alpha;
dEpsU    = -dEpsN;      % unloading strain increment
dEpsUU   =1*dEpsN;      % reloading strain increment
done     = 0;          % allow unloading once

noSteps1  = (epsilonr-epsilonF)/dEpsN;  % # of increments from peak down to unloading point
noSteps   = epsilonr/dEpsN;
noStepsF  = epsilonF/dEpsE; % number of elastic increments

% stored variables for plotting stress-strain...
stressA  = [0];
strainA  = [0];
jumpA    = [];
tracA    = [];

done1=0;
%% load increment loop
for i=1:noStepsF+noSteps1+noSteps+noSteps+noSteps
  %disp (sprintf('%s %i', 'Load increment: ',i) );
  
  if (i<=noStepsF) 
    dEps=dEpsE;
  elseif (done1==0)
    dEps = dEpsN;
  end
  
  if (i>noStepsF+noSteps1) && (done==0) 
    disp('unloading')
    i
    dEps = dEpsU; 
    done = 1;done1=1;
  end
  
  history.H       = H;
  history.cracked = status;
  history.loading = loading0;
  history.normal  = normal;
  history.kappa0  = kappa0;
  history.jump0   = jump0;
  history.ai      = ai;
  history.sigma0  = stress0;
  history.a       = a0     ;
  mat.ft          = mat.ft;
  
  load.eps0       = strain0;
  load.dEps       = [dEps;0;0];
  eps             = strain0+[dEps;0;0];
  
  %out             = updateTwoScaleCohesive(history,load,option);
  out             = updateTwoScaleCohesiveRK(history,load,option);
  
  ai              = reshape(out.K,1,4);
  sigma           = out.sigma;
  kappa           = out.kappa;
  jump            = out.u;
  loading         = out.loading;
  a               = reshape(out.tangent,1,9);
  
  % check failure
  sigmai = getPrincipalStress(sigma,mat.nu,mat.stressState);
  yield  = sigmai(1) - mat.ft;
  ft     = mat.ft;
  if (yield>=0)
    status        = 1;
    ai            = [-ft^2/mat.Gf 0 0 mat.ks];
    disp('cracking')
  end
  
    
  if jump(1)<0
    disp('negative jump!!!'); 
  end
  
  % store values for plotting
  stressA  = [stressA;sigma(1)];
  strainA  = [strainA;eps(1)];
  
  jumpA   = [jumpA;jump(1) jump(2)];
  tracA   = [tracA;out.trac(1) out.trac(2)];
  
  % update old quantities
  
  stress0  = sigma;
  strain0  = eps;
  kappa0   = kappa;
  jump0    = jump;
  loading0 = loading;
  ai0      = ai;
  a0       = a;  
  
%   if (i==noStepsF+noSteps1+noSteps+10)
%     dEps=dEpsUU; 
%     done1=1;
%     disp('load reversal')
%   end

end

%save('implicit.mat','stressA','strainA');
%h=load('implicit.mat');

% plot

figure(1)
set(gca,'FontSize',14)
hold on
plot(strainA,stressA,'red-*','LineWidth',1.9);
%plot(h.strainA,h.stressA,'black--','LineWidth',1.2);
%plot(h5.strainA,h5.stressA,'r*--','LineWidth',1.2);
%plot(h10.strainA,h10.stressA,'b*-','LineWidth',1.2);
%plot(h20.strainA,h20.stressA,'m*-','LineWidth',1.2);
%plot(strain0,stress,'black-','LineWidth',1.6);
%plot(jump,trac,'r-','LineWidth',1.6);
xlabel('strain')
ylabel('stress')
legend('RK','implicit','explicit,n=20')
axis([0.0 1.01 -0.0 1])

figure(2)
set(gca,'FontSize',14)
hold on
plot(jumpA(:,1),tracA(:,1),'b-','LineWidth',1.6);
%plot(strain0,stress,'black-','LineWidth',1.6);
%plot(jump,trac,'r-','LineWidth',1.6);
xlabel('jump')
ylabel('traction')
% legend('\alpha=1','\alpha=100','\alpha=1000','\alpha=1000,explicit')
%axis([0 8 0 1])

%save('tryH4Exp1000.mat','stress','strain','jump','trac');

% h1=load('tryH1Exp.mat');
% h2=load('tryH2Exp.mat');
% h4=load('tryH4Imp1.mat');
% h8=load('tryH10Exp.mat');

% h1=load('tryH4Imp1.mat');
% h100=load('tryH4Imp100.mat');
% h1000=load('tryH4Imp1000.mat');
% h1000Exp=load('tryH4Exp1000.mat');
%
% figure
% set(gca,'FontSize',14)
% hold on
% plot(h1.strain,h1.stress,'b-','LineWidth',1.6);
% plot(h100.strain,h100.stress,'r-','LineWidth',1.6);
% plot(h1000.strain,h1000.stress,'black-','LineWidth',1.6);
% plot(h1000Exp.strain,h1000Exp.stress,'cy-','LineWidth',1.6);
% %plot(strain0,stress,'black-','LineWidth',1.6);
% %plot(jump,trac,'r-','LineWidth',1.6);
% xlabel('strain')
% ylabel('stress')
% legend('\alpha=1','\alpha=100','\alpha=1000','\alpha=1000,explicit')
% axis([0 3 0 3])

% figure
% set(gca,'FontSize',14)
% hold on
% plot(h1.strain,h1.stress,'b-','LineWidth',1.6);
% plot(h2.strain,h2.stress,'r-','LineWidth',1.6);
% plot(h4.strain,h4.stress,'black-','LineWidth',1.6);
% plot(h8.strain,h8.stress,'black--','LineWidth',1.6);
% %plot(strain0,stress,'black-','LineWidth',1.6);
% %plot(jump,trac,'r-','LineWidth',1.6);
% xlabel('strain')
% ylabel('stress')
% legend('H=1','H=2','H=4','H=11')
%
% %%
%
% figure
% set(gca,'FontSize',14)
% hold on
% plot(h1.jump,h1.trac,'b.-','LineWidth',1.6);
% plot(h2.jump,h2.trac,'r--','LineWidth',1.6);
% plot(h4.jump,h4.trac,'cy*','LineWidth',0.3);
% plot(h8.jump,h8.trac,'black*','LineWidth',0.3);
% %plot(strain0,stress,'black-','LineWidth',1.6);
% %plot(jump,trac,'r-','LineWidth',1.6);
% xlabel('jump')
% ylabel('traction')
% legend('H=1','H=2','H=4','H=11')


