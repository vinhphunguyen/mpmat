function out=updateTwoScaleCohesiveRK(history,load,option)
% Stress update of the two-scale cohesive model.
% Explicit/Implicit stress update with Runge-Kutta for the explicit step.
% Bulk: elastic
% Crack: initially rigid exponential cohesive law.
% VP Nguyen, nvinhphu@gmail.com
% University of Adelaide, Australia, August, 2014.
%
% Status:
%
% 9/9/2014: add unloading with secant stiffness (Fig.6 in the
% 2scale-cohesive.pdf)
% 11/9/2014: add implicit step corrects unloading branch

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
%a       = reshape(history.a, 3,3);
u0      = history.jump0;   % jump in local coords of previous converged step
kappa0  = history.kappa0;  % history variable (maximum crack opening ever reached)
sigma   = history.sigma0;  % stress of previous converged step
loading0 = history.loading;% loading or not of previous step (increment)

sigma0=sigma;

% get loading
dEps    = load.dEps;      % strain increment
eps0    = load.eps0;      % strain at old step/increment

% implicit or explicit stress update
implicit  = option.implicit;
tolerance = option.tolerance;
tangent   = option.tangent; % need tangent stiffness or not
iterMax   = option.iterMax;
%stepCount = option.stepCount;

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
  
  rhs     = n'*(ao*dEps);
  
  C       = R'*K*R + Ao;
  DeltaU1 = inv(C)*rhs;           % jump increment
  u2      = u0 + R*DeltaU1*0.5;   % total jump in local coords.
  tsl     = updateExponentialTSL(kappa0,loading0,ft,Gf,ks,kp,u2);

  C       = R'*tsl.K*R + Ao;
  DeltaU2 = inv(C)*rhs;           % jump increment
  u3      = u0 + R*DeltaU2*0.5;   % total jump in local coords.
  tsl     = updateExponentialTSL(kappa0,loading0,ft,Gf,ks,kp,u3);
  C       = R'*tsl.K*R + Ao;
  DeltaU3 = inv(C)*rhs;           % jump increment
  u4      = u0 + R*DeltaU3;       % total jump in local coords.
  tsl     = updateExponentialTSL(kappa0,loading0,ft,Gf,ks,kp,u4);
  C       = R'*tsl.K*R + Ao;
  DeltaU4 = inv(C)*rhs;           % jump increment
  
  DeltaU  = (1/6)*( DeltaU1 + 2*DeltaU2 + 2*DeltaU3 + DeltaU4 );
  u       = u0 + R*DeltaU;
  tsl     = updateExponentialTSL(kappa0,loading0,ft,Gf,ks,kp,u);
  
    
%   if (tsl.comp)
%     disp('compression')
%     out.K       = tsl.K;
%     out.u       = u0;
%     out.sigma   = ao* (eps0 + dEps);%ao* (eps0 + dEps)
%     out.kappa   = tsl.kappa;
%     out.trac    = tsl.t;
%     out.loading = tsl.loading;
%     out.damage  = tsl.damage;
%     out.cracked = 1;
%     out.tangent = ao ;
%     return;
%   end
  
  % update total stress
  DeltaS  = ao*(dEps-(1/H)*n*DeltaU);
  sigma   = sigma + DeltaS;
  
  % the following is needed for resolving the load increment that 
  % initiated some cracks, dEps=0 and the iterative implicit
  % will not converge then.
  
  if norm(dEps) < 1e-16, implicit=0;end
  
  if (implicit)
    iterCount = 1;%sigma(1),tsl.t(1)
    residual  = n'*sigma-R'*tsl.t;    
    err       = norm(residual)/norm(n'*sigma);    
    %disp (sprintf('   %s %i %s %5.4e ', 'Implicit Stress, Iter',iterCount, ':', err) );
    while(err > tolerance)&&(iterCount<iterMax)
      C      = R'*K*R + Ao;
      deltaU = inv(C)*residual;
      u      = u + R*deltaU;      
      % update cohesive law      
      tsl      = updateExponentialTSL(kappa0,loading0,ft,Gf,ks,kp,u);
      deltaS   = -ao*(1/H)*n*deltaU;
      sigma    = sigma + deltaS;      
      K        = tsl.K;
      residual = n'*sigma-R'*tsl.t;
      err      = norm(residual)/norm(n'*sigma);      
      iterCount = iterCount + 1;
      %disp (sprintf('   %s %i %s %5.4e ', 'Implicit Stress, Iter',iterCount, ':', err) );
    end
    
    if (iterCount==iterMax) && (err > tolerance)
      tsl.comp, u0, u, tsl.t(1),sigma0(1),sigma(1)
      str=sprintf('   %s%d%s%d%s','Implicit stress update did not converge in ',iterCount, ' steps ',err,'\n');
      disp(str);
    end
    
  end
  
  sigma(sigma>0)=0; % set tensile stress to zero
  
  %write to output
  out.K       = tsl.K;
  out.u       = u;
  out.sigma   = sigma;
  out.kappa   = tsl.kappa;
  out.trac    = tsl.t;
  out.loading = tsl.loading;
  out.damage  = 1;
  out.cracked = 1;
  C           = R'*tsl.K*R + Ao;
  tangent     = ao - (1/H)*ao*(n*inv(C))*(n'*ao);
%   if (tsl.comp)
%     tangent=ao;
%     disp('compression, final')
%     %dEps,eps0
%   end
  out.tangent = tangent;
end


















