function [phi] = mlsConstantBasis1D(pt,index,node,di,form)
% Compute the MLS shape function at point pt for all nodes within the
% support of this point pt.
% Basis used is linear basis pT = [1]
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 9 October 2015.

% --------------------------------------
%      compute the moment matrix A
% --------------------------------------
ni   = length(index);
A    = 0;
phi  = zeros(ni,1);
w    = zeros(ni,1);

for m = 1 : ni
  idx  = index(m);
  xi   = node(idx);
  wi   = computeCircleSpline(pt,xi,di(idx),form);  
  A    = A    + wi;
  w(m) = wi ;
end

for m = 1 : ni
  xi  = node(index(m));  
  phi(m) = w(m)/A;
end
