function [phi] = mlsLinearBasis2D(pt,index,node,di,form)
%#codegen
% Compute the MLS shape function at point pt for all nodes within the
% support of this point pt.
% Basis used is linear basis pT = [1 x y]
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 9 October 2015.

% --------------------------------------
%      compute the moment matrix A
% --------------------------------------
nn   = length(index);
A    = zeros(3,3);
w    = zeros(nn,1);
phi  = zeros(nn,1);

for m = 1:nn
  idm  = index(m);
  xi   = node(idm,1);
  yi   = node(idm,2);
  wi   = computeCircleSpline(pt,[xi yi],di(idm),form);
%   wi   = 0;
%   r    = norm( pt - [xi yi] )/di(idm);
%   if ( r <= 1. )
%     wi = 1 - 6*r^2 + 8*r^3 -3*r^4;
%   end
  pTp  = [1 xi yi;xi xi^2 xi*yi;yi xi*yi yi^2] ;
  A    = A    + wi*pTp ;
  w(m) = wi;
end

p  = [1; pt(1); pt(2)];
c = A\p;

for m = 1:nn
  xi  = node(index(m),:);
  piT = [1; xi(1); xi(2)];
  phi(m) = c'*piT*w(m) ;
end

%codegen -o mlsLinearBasis2D_mex mlsLinearBasis2D -args {coder.typeof(double(0), [1 2]),coder.typeof([1 2],[Inf 1]),coder.typeof(double(0), [Inf 2]),coder.typeof(double(0), [Inf 1]),coder.typeof('aa')}
