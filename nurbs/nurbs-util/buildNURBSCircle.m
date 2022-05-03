function solid = buildNURBSCircle(r,center)

% knots
uKnot = [0 0 0 1 1 1];
vKnot = [0 0 0 1 1 1];

% control points for r=0.5 and (0,0);
controlPts          = zeros(4,3,3);

controlPts(1:2,1,1) = [-sqrt(2)/4 sqrt(2)/4];
controlPts(1:2,2,1) = [-sqrt(2)/2 0];
controlPts(1:2,3,1) = [-sqrt(2)/4 -sqrt(2)/4];

controlPts(1:2,1,2) = [0 sqrt(2)/2];
controlPts(1:2,2,2) = [0;0];
controlPts(1:2,3,2) = [0 -sqrt(2)/2];

controlPts(1:2,1,3) = [sqrt(2)/4 sqrt(2)/4];
controlPts(1:2,2,3) = [sqrt(2)/2 0];
controlPts(1:2,3,3) = [sqrt(2)/4 -sqrt(2)/4];

controlPts(1:2,:,:) = r/0.5*controlPts(1:2,:,:);

% weights
controlPts(4,:,:)   = 1;

fac = sqrt(2)/2;

controlPts(4,2,1) = fac;
controlPts(4,1,2) = fac;
controlPts(4,3,2) = fac;
controlPts(4,2,3) = fac;

controlPts(1:2,2,1) = fac*controlPts(1:2,2,1);
controlPts(1:2,1,2) = fac*controlPts(1:2,1,2);
controlPts(1:2,3,2) = fac*controlPts(1:2,3,2);
controlPts(1:2,2,3) = fac*controlPts(1:2,2,3);

% build NURBS object
solid = nrbmak(controlPts,{uKnot vKnot});
solid = nrbtform(solid,vectrans(center));
