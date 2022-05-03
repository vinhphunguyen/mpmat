function [bx,by] = vortexBodyForces(x,y,time,mu,rho0,G,Ri,Ro)
%
% MMS body force for vortex problem.
%
% x,y: MUST BE initial coords!!!
%

T     = 1;
r     = sqrt(x*x+y*y);
theta = atan2(y,x);
R     = (Ri+Ro)/2;
s     = (r - R)/(Ri-Ro);
h     = 1 - 8*((r - R)/(Ri-Ro))^2 +16*((r - R)/(Ri-Ro))^4;
hp    = - 16*(r - R)/(Ri-Ro)^2 + 16*4*(r - R)^3/(Ri-Ro)^4;
hpp   = - 16/(Ri-Ro)^2 + 16*4*3*(r - R)^2/(Ri-Ro)^4;
g     = G * sin(pi*time/T);
gp    = G*pi/T*cos(pi*time/T);
gpp   = -pi*pi/(T*T)*g;
alpha = g*h;
mdr   = mu/rho0;

br    = ( mdr*(3*g*hp+r*g*hpp) - r*gpp*h)*sin(alpha) + (mdr*r*(g*hp)^2 - r*(gp*h)^2)*cos(alpha);
bt    = (-mdr*(3*g*hp+r*g*hpp) + r*gpp*h)*cos(alpha) + (mdr*r*(g*hp)^2 + r*(gp*h)^2)*sin(alpha);

bx    = br*cos(theta) - bt*sin(theta);
by    = br*sin(theta) + bt*cos(theta);



