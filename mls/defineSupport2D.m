function [index] = defineSupport2D(node,x,di)
% find nodes in neighbouring of point x
% Inpputs:
%  node : numnode x 2, nodal coordinates
%  x    : 1 x 2, coordinate of point 
%  di   : 1 x numnode, size of support of nodes
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 9 October 2015.

numnode = size(node,1) ;
dif     = node - [ones(numnode,1)*x(1) ones(numnode,1)*x(2)];
r       = sqrt(sum(dif.^2,2));
index   = find(r - di <= 0.00001);
