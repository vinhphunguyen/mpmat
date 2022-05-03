function [ra,rb] = vortexRadii (Xa,Ya,Xb,Yb,t,A)
% compute the updated inner/outer radii for points (Xa,Ya) and (Xb,Yb).
%

Ra = sqrt(Xa^2+Ya^2);
Rb = sqrt(Xb^2+Yb^2);

alphaA = A*sin(pi*t)*(1-32*(Ra-1)^2+256*(Ra-1)^4);
alphaB = A*sin(pi*t)*(1-32*(Rb-1)^2+256*(Rb-1)^4);

Qa = [cos(alphaA) -sin(alphaA) 0;...
      sin(alphaA)  cos(alphaA) 0;...
      0            0          1];
    
Qb = [cos(alphaB) -sin(alphaB) 0;...
      sin(alphaB)  cos(alphaB) 0;...
      0            0          1];    
    
    
xa = Qa*[Xa;Ya;0];    
xb = Qb*[Xb;Yb;0];

ra = sqrt(xa(1)^2+xa(2)^2);
rb = sqrt(xb(1)^2+xb(2)^2);

