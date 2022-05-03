function [phi] = mlsLinearBasis1D(pt,index,node,di,form)
% Compute the MLS shape function at point pt for all nodes within the
% support of this point pt.
% Basis used is linear basis pT = [1 x]
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 9 October 2015.

% --------------------------------------
%      compute the moment matrix A
% --------------------------------------
ni   = length(index);
A    = zeros(2,2);
phi  = zeros(ni,1);
w    = zeros(ni,1);

for m = 1 : ni
  idx  = index(m);
  xi   = node(idx);
  wi   = computeCircleSpline(pt,xi,di(idx),form);
  pTp  = [1 xi]'*[1 xi] ;
  A    = A    + wi*pTp ;
  w(m) = wi ;
end

p  = [1; pt];
% --------------------------------------
%         compute  matrix c(x)
% --------------------------------------
% A(x)c(x)   = p(x)
% Backward substitutions, two times for c(x)

c = A\p;

for m = 1 : ni
  xi  = node(index(m));
  piT = [1;xi];
  phi(m) = c'* piT*w(m);
end
