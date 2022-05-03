% An example of the two-scale damage model.
% Bulk: elastic
% Crack: isotropic damage model.
%
% Note: normal to the localisation band is known a priori.
%
% VP Nguyen
% University of Adelaide, Australia, August, 2014.

addpath ../fem-functions/

clear all;
clc;

E  = 5;
nu = 0.;
ft = 1;
Gf = 1;
ki  = ft/E;       % damage threshold
kc  = 4;          % critical damage (used in linear softening law) 
alpha=0.999;      % constant in damage evolution law
beta = 1;        % constant in damage evolution law

h   = 0.5;       % internal length scale
H   = 2;
eta = h/H;

cracked = 0;

De = elasticityMatrix(E,nu,'PLANE_STRAIN');
n  = [1 0; 0 0; 0 1];
Ao = n'*De*n;
ao = De;

nstep = 100;

epsilon  = 0; sigma=0; epsi = [0;0;0]; epsilon0=0;
epsilonf = 1;
dEps     = 5e-2; %1e-4, 1e-3, 1e-2: ok; 5e-2: no ok
kappa0   = 0;

% stress-strain of the whole  
stress  = [];
strain  = [];
strain0 = [];

% stress-strain of localisation band
stressi   = [];
straini   = [];

while (epsilon < epsilonf)
  if (cracked==0)
    eps = [epsilon;0;0];
    eqvS = sqrt((1./E)*eps'*De*eps);
    f    = eqvS - kappa0;
    kappa = kappa0;
    if ( f > 0.) 
      kappa = eqvS;
      loading=1;
    else
      loading=0;
    end
    damage = 0.;
    if ( kappa >= ki )
      damage = 1. - (ki/kappa)*(1.-alpha+alpha*exp(-beta*(kappa-ki)));
    end
    
    sigma = (1.-damage)*De*eps;
    
    % compute ai
    if (loading) && (damage)
    a        = ki / (kappa * kappa);
    b        = a * kappa;
    c        = alpha * exp(-beta * (kappa-ki) );
    dddk     = a * (1.0 - alpha+ c) + beta * b * c;
    
    deqvdeps = (1/E/eqvS)*De*eps;
    ai       = (1-damage)*De - loading*dddk*De*eps*deqvdeps'
    
    ao       = (1-damage)*De;
    Ait      = n'*ai*n;
    dd=det(Ait)
    det(ai)
    % check localisation occurence
    if (dd<=0)
      cracked = 1;
      epsilon, damage
      epsi(1)  = epsilon;     
    end
    end
  else
    Ao   = n'*ao*n;
    Ai   = n'*ai*n;
    C        = (eta/h)*Ao + ((1-eta)/h)*Ai;
    djumpu   = inv(C)*n'*ao*[dEps;0;0];
    if (djumpu(1)<0)
      (eta/h)*Ao
      ((1-eta)/h)*Ai
      error('negative jump')
    end
    depsi    = (1/h)*n*djumpu;
    % update strain at the localisation band
    epsi = epsi + depsi;
    % update localisation band
    
    eqvS = sqrt((1./E)*epsi'*De*epsi);
    f    = eqvS - kappa0;
    kappa = kappa0;
    if ( f >= 0.) kappa = eqvS; end
    damage = 0.;
    if ( kappa >= ki )
      damage = 1. - (ki/kappa)*(1.-alpha+alpha*exp(-beta*(kappa-ki)));
    end        
    sigmai = (1-damage)*De*epsi;
    % compute ai
    a        = ki / (kappa * kappa);
    b        = a * kappa;
    c        = alpha * exp(-beta * (kappa-ki) );
    dddk     = a * (1.0 - alpha+ c) + beta * b * c;
    
    deqvdeps = (1/E/eqvS)*De*epsi;
    ai       = (1-damage)*De - dddk*De*epsi*deqvdeps';
    
    % compute stress increment at outside (also the total)
    depso    = [dEps;0;0] - eta*depsi;
    dsigma   = (1/(1-eta))*ao*depso;
    sigma    = sigma + dsigma(1);
    
    stressi = [stressi;sigmai(1)];
    straini = [straini;epsi(1)];
  end
  
  stress = [stress;sigma(1)];
  strain = [strain;epsilon];
  strain0 = [strain0;epsilon0];
  
  epsilon = epsilon + dEps;
  kappa0 = kappa;
end

save('tryH2.mat','stress','strain','stressi','straini');

h1=load('tryH1.mat');
h2=load('tryH2.mat');
h4=load('tryH4.mat');

figure
set(gca,'FontSize',14)
hold on
plot(h1.strain,h1.stress,'b-','LineWidth',1.6);
plot(h2.strain,h2.stress,'r-','LineWidth',1.6);
%plot(h4.strain,h4.stress,'black-','LineWidth',1.6);
% plot(h4.strain,h4.stress,'black-','LineWidth',1.6);
% %plot(strain0,stress,'black-','LineWidth',1.6);
% %plot(jump,trac,'r-','LineWidth',1.6);
xlabel('averaged strain')
ylabel('averaged stress')
legend('H=1','H=2')


figure
set(gca,'FontSize',14)
hold on
plot(h1.straini,h1.stressi,'b.','LineWidth',1.6);
plot(h2.straini,h2.stressi,'r-','LineWidth',1.6);
%plot(h4.straini,h4.stressi,'black-','LineWidth',1.6);
% plot(h4.strain,h4.stress,'black-','LineWidth',1.6);
% %plot(strain0,stress,'black-','LineWidth',1.6);
% %plot(jump,trac,'r-','LineWidth',1.6);
xlabel('localised strain')
ylabel('localised stress')
legend('H=1','H=2')


