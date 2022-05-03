function [w,dwdr] = cubicSpline(r)
% Compute cubic spline function

if (r <= 0.5)
    w    = 2/3 - 4*r*r + 4*r*r*r ;
    dwdr = -8*r + 12*r*r ;
elseif (r > 0.5) && (r <= 1.0)
    w    = 4/3 - 4*r + 4*r*r - (4/3)*r*r*r ;
    dwdr = -4  + 8*r - 4*r*r ;
else
    w    = 0.0 ;
    dwdr = 0.0 ;
end
