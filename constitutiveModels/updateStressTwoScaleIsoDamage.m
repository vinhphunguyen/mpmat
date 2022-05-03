function [out] = updateStressTwoScaleIsoDamage(dstrain,mat,hist,straini)
%
% update stress for isotropic elastic-damage model in a two-scale
% constitutive model of GD Nguyen.
% Inputs:
%  strain:  total strain of a particle
%  matProp: a structure containing material properties
%  kappa0: internal variables, maximum eqv. strain ever reached.
% Outputs: updated stress
% VP Nguyen
% Saigon, Vietnam, July 2014.

E   = mat.young;
De  = mat.De;
ki  = mat.ki;
al  = mat.alpha;
bt  = mat.beta;

kappa0     = hist.kappa;
localised  = hist.localised;
strain0    = hist.strain;

if (localised==0)
 
  out = updateStressIsoDamage(strain,mat,kappa0,loading0);
  
else
  h  = matProp.bandWidth;
  H  = matProp.elemWidth;
  eta=h/H;
  
  ao = hist.ao;
  strain = straini;
  
  eqv = sqrt((1./E)*strain'*De*strain);
  f   = eqv - kappa0;
  
  kappa = kappa0;
  
  if ( f >= 0.)
    kappa = eqv;
  end
  damage = 0.;
  if ( kappa >= ki )
    damage = 1. - (ki/kappa)*(1.-al+al*exp(-bt*(kappa-ki)));
  end
 
  a        = Ki / (K * K);
  b        = a * K;
  c        = alpha * exp(-beta * (kappa-ki) );
  dddk     = a * (1.0 - al+ c) + bt * b * c;
  
  eqv      = sqrt((1./E)*strain'*De*strain);
  deqvdeps = (1/E/eqv)*De*strain;
  ai       = (1-damage)*De - dddk*De*strain*deqvdeps';
  
  nx       = 1; ny = 0; % hard coded
  n        = zeros(3,2);
  n(1,1)   = nx;
  n(2,2)   = ny;
  n(3,1)   = ny;
  n(3,2)   = nx;
  Ao       = n'*ao*n;
  Ai       = n'*ai*n;
  C        = (eta/h)*Ao + (1-eta/h)*Ai;
  
  djumpu   = inv(C)*n'*ao*dstrain;
  depsi    = (1/h)*n*djumpu;
  
  straini = straini + depsi;
  
  depso    = dstrain - eta*depsi;
  dsigma   = (1/(1-eta))*ao*depso;
  stress   = stress + dsigma;
end

out.stress = stress;
out.localised = localised;
out.kappa     = kappa;



