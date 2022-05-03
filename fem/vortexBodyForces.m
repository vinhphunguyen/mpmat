function [bx,by] = vortexBodyForces(x,y,time,mu,rho0,G,Ri,Ro)
%
% MMS body force for vortex problem.
%
% x,y: MUST BE initial coords!!!
%
% r     = sqrt(x*x+y*y);
% theta = atan2(y,x);
% p1    = 4096*r*(15-47*r+48*r^2-16*r^3)^2*mu*(sin(2*pi*t))^4;
% p1    = p1/rho0;
% p2    = pi*pi*r*(15-32*r+16*r^2)^4*(sin(2*pi*t))^2;
% p3    = -16*(-45+188*r-240*r^2+96*r^3);
% p4    =     (-45+188*r-240*r^2+96*r^3);
% p5    =     (15-32*r+16*r^2)^2;
% br    = p1 - p2;
% bt    = (1/rho0)*( 2*mu*p3 + 2*cos(2*pi*t)*(16*mu*p4+pi^2*r*rho0*p5) );
% % (15-32*R+16*R^2)^2 equlas (1-32*(R-1)^2+256*(R-1)^4)!!!
% %alpha = A*0.5*(1-cos(2*pi*t))*(15-32*r+16*r^2)^2;
% alpha = A*0.5*(1-cos(2*pi*t))*(1-32*(r-1)^2+256*(r-1)^4);
% % Uintah implementation (requires initial u/v)
% % br=-r*A^2*(cos(pi*t))^2*pi^2*(4*r-3)^4*(4*r-5)^4+...
% %  (2/rho0)*2048*(4*r-5)^2*(r-1)^2*(4*r-3)^2*A^2*(sin(pi*t))^2*r*mu;
% % bt=-A*sin(pi*t)*(r*pi^2*(4*r-3)^2*(4*r-5)^2+...
% %   (1/rho0)*64*mu*(96*r^3-240*r^2+188*r-45));
% % alpha = A*sin(pi*t)*(1-32*(r-1)^2+256*(r-1)^4);
% thetan = theta + alpha;
% bx    = br*cos(thetan) - bt*sin(thetan);
% by    = br*sin(thetan) + bt*cos(thetan);

% From Alban
PI    = pi;
T     = 1;
r     = sqrt(x*x+y*y);
theta = atan2(y,x);
R     = (Ri+Ro)/2;
s     = (r - R)/(Ri-Ro);
h     = 1 - 8*((r - R)/(Ri-Ro))^2 +16*((r - R)/(Ri-Ro))^4;
hp    = - 16*(r - R)/(Ri-Ro)^2 + 16*4*(r - R)^3/(Ri-Ro)^4;
hpp   = - 16/(Ri-Ro)^2 + 16*4*3*(r - R)^2/(Ri-Ro)^4;
g     = G * sin(PI*time/T);
gp    = G*PI/T*cos(PI*time/T);
gpp   = -PI*PI/(T*T)*g;
alpha = g*h;
mdr   = mu/rho0;

br    = ( mdr*(3*g*hp+r*g*hpp) - r*gpp*h)*sin(alpha) + (mdr*r*(g*hp)^2 - r*(gp*h)^2)*cos(alpha);
bt    = (-mdr*(3*g*hp+r*g*hpp) + r*gpp*h)*cos(alpha) + (mdr*r*(g*hp)^2 + r*(gp*h)^2)*sin(alpha);

bx    = br*cos(theta) - bt*sin(theta);
by    = br*sin(theta) + bt*cos(theta);



