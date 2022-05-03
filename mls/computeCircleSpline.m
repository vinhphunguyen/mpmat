function [w] = computeCircleSpline(x,xI,d,form)
% Compute cubic and quartic spline function
% Inputs:
% x (1x2)  : coordinate of point at which w is to be evaluated
% xI (1x2) : coord of node I
% d        : size of the support

r = norm( x - xI )/d;

switch form
  case 'cubic_spline' 
     [w,~] = cubicSpline(r);
  case 'quartic_spline'
     [w,~] = quarticSpline(r);
  otherwise 
     error('Grr. Unknown functional form');
end


