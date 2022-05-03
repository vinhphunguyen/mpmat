function [w,dwdr] = expSpline(r,alpha)
% Compute exponential spline function

if (r <= 1.)
    w    = exp(-(r/alpha)^2) ;
    dwdr = -8*r + 12*r*r ;
else
    w    = 0.0 ;
    dwdr = 0.0 ;
end
