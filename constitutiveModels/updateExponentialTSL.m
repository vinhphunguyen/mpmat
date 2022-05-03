function out=updateExponentialTSL(kappa0,loading0,ft,Gf,ks,kp,u)
% Stress update of the exponential cohesive model of G.N. Wells and L.J. Sluys 
% (IJNME,2001), New finite elements to cohesive crack modeling.
% VP Nguyen, nvinhphu@gmail.com/ phu.nguyen@adelaide.edu.au
% University of Adelaide, Australia, August, 2014.

comp    = 0;
uc      = u;
if (uc(1)<0)
  uu=uc(1);uc(1)=0; comp=1;
end

%kappa = norm(uc);
kappa = uc(1);

f = kappa - kappa0;                     % loading function

if     (abs(f) < 1e-14 )
  loading = loading0;
elseif ( f < 0 )                           % unloading
  loading = 0;  
else                                       % loading
  loading = 1;
end

if     (loading)                           % loading
  tn    = ft*exp(-(ft/Gf)*kappa);
  ts    = ks*u(2);
  K11   = -ft^2/Gf*exp(-(ft/Gf)*kappa);
  K22   = ks;
  damage=1-tn/ft;
elseif (abs(kappa0)<1e-8)                  % unloading
  %disp('dd')
  K11  = -ft^2/Gf;
  tn   = ft;
  ts   = ks*u(2);
  K22  = ks;
  damage=0;
else                                       % secant unloading
  t0   = ft*exp(-(ft/Gf)*kappa0);
  K11  = t0/kappa0;
  tn   = K11*u(1);%K11*kappa;
  ts   = ks*u(2);
  K22  = ks;
  damage=1-t0/ft;                          % damage frozen at previous step
end

if (f<0) kappa   = kappa0; end

%penalty to prevent crack penetration

if comp
  %disp('com')
  tn = kp*u(1);
  K11= kp;
end

if (kappa0==0) && (kappa==0)
  %disp('fff'),comp
  K11  = -ft^2/Gf;
  tn   = ft;
  ts   = ks*u(2);
  K22  = ks;
  damage=0;
end

out.loading = loading;
out.t       = [tn; ts];
out.K       = [K11 0;0 K22];
out.damage  = damage;
out.kappa   = kappa;
out.comp    = comp;

