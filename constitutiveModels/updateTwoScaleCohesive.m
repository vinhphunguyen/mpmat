function out=updateTwoScaleCohesive(history,load,option)
% Stress update of the two-scale cohesive model.
% Explicit/Implicit stress update.
% Bulk: elastic
% Crack: initially rigid exponential cohesive law.
% VP Nguyen, nvinhphu@gmail.com
% University of Adelaide, Australia, August, 2014.
%
% Status:
%
% 9/9/2014: add unloading with secant stiffness (Fig.6 in the
% 2scale-cohesive.pdf)

global mat

% get material properties
ft      = mat.ft;
Gf      = mat.Gf;
ao      = mat.ao;
kp      = mat.penalty;
ks      = mat.ks;

% get history variables
H       = history.H;
cracked = history.cracked;
normal  = history.normal; % vector
K       = reshape(history.ai,2,2);
a       = reshape(history.a, 3,3);
u0      = history.jump0;   % jump in local coords of previous converged step
kappa0  = history.kappa0;  % history variable (maximum crack opening ever reached)
sigma   = history.sigma0;  % stress of previous converged step
loading0 = history.loading;% loading or not of previous step (increment)

% get loading
dEps    = load.dEps;      % strain increment
eps0    = load.eps0;      % strain at old step/increment

% implicit or explicit stress update
implicit  = option.implicit;
tolerance = option.tolerance;
tangent   = option.tangent; % need tangent stiffness or not
iterMax   = option.iterMax;
stepCount = option.stepCount;

% compute normal in matrix form, 1/H*Ao and rotation matrix R
n       = [normal(1) 0;0 normal(2); normal(2) normal(1)];
Ao      = (1/H)*n'*ao*n;
R       = [normal(1) normal(2);-normal(2) normal(1)];

%% elastic stage
if (cracked==0)
  sigma       = ao * (eps0 + dEps);
  out.sigma   = sigma;
  out.tangent = ao;
  out.cracked = 0;
  out.K       = [0 0; 0 0];
  out.u       = [0; 0];
  out.kappa   = 0;
  out.trac    = [0;0];
  out.loading = 0;
  out.damage  = 0;
  out.cracked = 0;
  %% cracked stage, two scale model
else
  u00=u0;u=u0;
  %1. EXPLICIT STEP
  dEpsi  = (dEps/stepCount);
  for i=1:stepCount
    C      = R'*K*R + Ao;
    DeltaU = inv(C)*n'*(ao*dEpsi);        % jump increment
    u      = u + R*DeltaU;                % total jump in local coords.
    %% update cohesive law
    tsl = updateExponentialTSL(kappa0,loading0,ft,Gf,ks,kp,u);
    % update total stress
    DeltaS = ao*(dEpsi-(1/H)*n*DeltaU);
    sigma  = sigma + DeltaS;
    %end
    K      = tsl.K; %K0 = K;
    tn     = tsl.t(1);
    ts     = tsl.t(2);%u0=u;
  end
  
  %2. IMPLICIT STEP or CORRECTION
  if (implicit)
    iterCount = 1;
    residual  = n'*sigma-R'*[tn;ts];
    err       = norm(residual)/norm(n'*sigma);
    %disp (sprintf('   %s %i %s %5.4e ', 'Implicit Stress, Iter',iterCount, ':', err) );
    while(err > tolerance)&&(iterCount<iterMax)
      C      = R'*K*R + Ao;
      deltaU = inv(C)*residual;
      u      = u + R*deltaU;      %u,u0
      % update cohesive law
      tsl    = updateExponentialTSL(kappa0,loading0,ft,Gf,ks,kp,u);
      %       if(tsl.loading==0)
      %         disp('negative strain increment');%tsl.t(1),u,tsl.loading,loading0,sigma
      %       end
      % update total stress
      deltaS = ao*(-(1/H)*n*deltaU);
      sigma  = sigma + deltaS;
      K        = tsl.K;
      residual = n'*sigma-R'*tsl.t;
      err      = norm(residual)/norm(n'*sigma);
      %err      = residual(1)/norm(n'*sigma);
      iterCount = iterCount + 1;
      %disp (sprintf('   %s %i %s %5.4e ', 'Implicit Stress, Iter',iterCount, ':', err) );
    end
    %str=sprintf('   %s%d%s%d','converged in ',iterCount, ' steps ',err);disp(str);
    if (iterCount==iterMax) && (err > tolerance)
      %tsl.comp
      str=sprintf('%s%d%s%d','Implicit stress update did not converge in ',iterCount, ' steps ',err);
      disp(str);
    end
  end
  
  out.K       = tsl.K;
  out.u       = u;
  out.sigma   = sigma;
  out.kappa   = tsl.kappa;
  out.trac    = tsl.t;
  out.loading = tsl.loading;
  out.damage  = tsl.damage;
  out.cracked = 1;
  C           = R'*K*R + Ao;
  tangent     = ao - (1/H)*ao*(n*inv(C))*(n'*ao);
  if (tsl.comp) tangent=ao;end
  out.tangent = tangent;
end


















