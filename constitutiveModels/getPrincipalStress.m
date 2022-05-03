function sigmai = getPrincipalStress(sigma,nu,stressState)
% Find principal stresses for a 2D stress state
% Both plane strain and plane stress is supported.
% VP Nguyen, nvinhphu@gmail.com
% University of Adelaide, Australia
% August 2014.

sigmaxx = sigma(1);
sigmayy = sigma(2);
sigmaxy = sigma(3);

d0      = (sigmaxx-sigmayy)^2+4.*sigmaxy^2;
d       = sqrt(d0);
temp    = 0.5*(sigmaxx+sigmayy);

sigmai(1) = temp + 0.5*d;
sigmai(2) = temp - 0.5*d;

if strcmp(stressState,'PLANE_STRAIN')
  sigmai(3) = nu*(sigmai(1)+sigmai(2));
end