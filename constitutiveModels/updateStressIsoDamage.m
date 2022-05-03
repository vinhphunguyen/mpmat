function out = updateStressIsoDamage(strain,mat,kappa0,loading0)
%
% update stress for isotropic elastic-damage model
% Inputs:
%  strain:  total strain of a particle
%  matProp: a structure containing material properties
%  kappa0: internal variables, maximum eqv. strain ever reached.
% Outputs: updated stress
% VP Nguyen
% Saigon, Vietnam, June 2014.

E   = mat.E;
De  = mat.ao;
ki  = mat.ki;
al  = mat.alpha;
bt  = mat.beta;

%eqv = sqrt((1./E)*strain'*De*strain);
eqv = computeMazarsEqvStrain(strain,mat);

f   = eqv - kappa0;


if     (abs(f) < 1e-14 )
  loading = loading0;
  kappa   = kappa0;
elseif ( f < 0 )                           % unloading
  loading = 0;
  kappa   = kappa0;
else                                       % loading
  loading = 1;
  kappa   = eqv;
end

damage = 0.;
if ( kappa > ki )
  %al,bt,ki,kappa,strain
  damage = 1. - (ki/kappa)*(1-al + al*exp(-bt*(kappa-ki)));
else
  loading=0;
end

if (abs(damage-1)<1e-9)
%  kappa,kappa0,eqv,strain
end

out.sigma  = (1-damage)*De*strain;
out.damage = damage;
out.loading=loading;
out.kappa  = kappa;
% 
if (loading)
  a        = ki / (kappa * kappa);
  b        = a * kappa;
  c        = al * exp(-bt * (kappa-ki) );
  dddk     = a * (1.0 - al+ c) + bt * b * c;
  
  %eqv      = sqrt((1./E)*strain'*De*strain);
  %deqvdeps = (1/E/eqv)*De*strain;
  
  deqvdeps = computeMazarsEqvStrainDer(strain,mat);
  
  out.tangent = (1-damage)*De - dddk*De*strain*deqvdeps';
  %out.tangent
else
  out.tangent = (1-damage)*De;
 end



