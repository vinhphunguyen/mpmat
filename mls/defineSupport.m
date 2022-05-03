function [index] = defineSupport(node,x,di)
% find nodes in neighbouring of point x
% Inpputs:
%  node : numnode x 1, nodal coordinates
%  x    : scalar, coordinate of point 
%  di   : 1 x numnode, size of support of nodes
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 9 October 2015.

numnode = size(node,1) ;
dif     = node - ones(numnode,1)*x;
r       = abs(dif);
index   = find(r - di <= 1e-16);
