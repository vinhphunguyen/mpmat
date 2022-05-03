function out=updateTwoScaleDamageRK(history,load,option)
% Stress update of the two-scale isotropic damage model.
% Explicit/Implicit stress update with Runge-Kutta for the explicit step.
% Bulk: elastic
% Crack: isotropic damage model
% VP Nguyen, nvinhphu@gmail.com
% University of Adelaide, Australia, September, 2014.
%
% Status:


global mat mesh

% get material properties
ft      = mat.ft;
Gf      = mat.Gf;
ao      = mat.ao;
h       = mat.h;


% get history variables
H       = history.H;
cracked = history.cracked;
normal  = history.normal; % vector
ai      = reshape(history.ai,3,3);
ao      = reshape(history.a ,3,3);
epsi0   = history.jump0;   % jump in local coords of previous converged step
kappa0  = history.kappa0;  % history variable (maximum crack opening ever reached)
sigma   = history.sigma0;  % stress of previous converged step
loading0 = history.loading;% loading or not of previous step (increment)

sigma0=sigma;

% get loading
dEps    = load.dEps;      % (averaged) strain increment
eps0    = load.eps0;      % (averaged) strain at old step/increment

% implicit or explicit stress update
implicit  = option.implicit;
tolerance = option.tolerance;
tangent   = option.tangent; % need tangent stiffness or not
iterMax   = option.iterMax;
stepCount = option.stepCount;

% compute normal in matrix form, 1/H*Ao and rotation matrix R
n       = [normal(1) 0;0 normal(2); normal(2) normal(1)];

%% elastic stage
%if (cracked==0)
  %eps0+dEps
  res = updateStressIsoDamage(eps0+dEps,mat,kappa0,loading0);
  %eps0,dEps
  out.sigma   = res.sigma;
  out.tangent = res.tangent;
  out.cracked = 0;
  out.K       = res.tangent;
  out.u       = eps0+dEps;
  out.kappa   = res.kappa;
  out.trac    = [0;0];
  out.loading = res.loading;
  out.damage  = res.damage;
  out.cracked = 0;
  
  %% cracked stage, two scale model
% else  
%   eta      = h/H;
%   Ao       = n'*ao*n;
%   
%   c1      = (1-eta)/h;
%   c2      = eta/h;
%   
%   %1. EXPLICIT STEP
%   
%   rhs     = n'*(ao*dEps);
%   C       = c1*(n'*ai*n) + c2*Ao;  
%   DeltaU1 = inv(C)*rhs;           % jump increment
%   depsi   = (1/h)*n*DeltaU1/2;
%   res     = updateStressIsoDamage(epsi0+depsi,mat,kappa0,loading0);
%             
%   C       = c1*(n'*res.tangent*n) + c2*Ao;
%   DeltaU2 = inv(C)*rhs;           % jump increment
%   depsi   = (1/h)*n*DeltaU2/2;
%   res     = updateStressIsoDamage(epsi0+depsi,mat,kappa0,loading0);
%   C       = c1*(n'*res.tangent*n) + c2*Ao;
%   DeltaU3 = inv(C)*rhs;           % jump increment
%   depsi   = (1/h)*n*DeltaU3;
%   res     = updateStressIsoDamage(epsi0+depsi,mat,kappa0,loading0);
%   C       = c1*(n'*res.tangent*n) + c2*Ao;   
%   DeltaU4 = inv(C)*rhs;           % jump increment
%   
%   DeltaU  = (1/6)*( DeltaU1 + 2*DeltaU2 + 2*DeltaU3 + DeltaU4 );
%   depsi   = (1/h)*n*DeltaU;
%   res     = updateStressIsoDamage(epsi0+depsi,mat,kappa0,loading0);
%   C       = c1*(n'*res.tangent*n) + c2*Ao;
%   % update total stress
%   depso    = dEps - eta*depsi;
%   dsigma   = (1/(1-eta))*ao*depso;
%   sigma    = sigma + dsigma;
%   
%   %write to output
%   out.K       = res.tangent;
%   out.u       = epsi0+depsi;
%   out.sigma   = sigma;
%   out.kappa   = res.kappa;
%   out.trac    = [0;0];
%   out.loading = res.loading;
%   out.damage  = res.damage;
%   out.cracked = 1; 
%   tangent     = (1/(1-eta))*(ao - (1/H)*ao*(n*inv(C))*(n'*ao));
%   out.tangent = tangent;
% end


















