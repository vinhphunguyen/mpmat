function [w,dwdr] = quarticSpline(r)
% Compute cubic spline function

if (r <= 1.0)
   w    = 1 - 6*r*r + 8*r*r*r - 3*r*r*r*r;
   dwdr = -12*r + 24*r*r - 12*r*r*r ;
else
   w    = 0.0 ;
   dwdr = 0.0 ;
end
