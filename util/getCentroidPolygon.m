function [centroid] = getCentroidPolygon (particles,p)

vx  = particles.node(particles.elem{p},1);
vy  = particles.node(particles.elem{p},2);
nv  = length(particles.elem{p});
vxS = vx([2:nv 1]);
vyS = vy([2:nv 1]); %Shifted vertices
temp       = vx.*vyS - vy.*vxS;
a          = 0.5*sum(temp);
centroid   = 1/(6*a)*[sum((vx+vxS).*temp),sum((vy+vyS).*temp)];
