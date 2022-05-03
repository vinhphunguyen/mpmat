%
%
%
%
% Vinh Phu Nguyen, March 2012
% nvinhphu@gmail.com

function mesh  = generateIGA2DMesh(nurbs)

uKnot          = nurbs.knots{1};
vKnot          = nurbs.knots{2};
p              = nurbs.order(1)-1;
q              = nurbs.order(2)-1;

uniqueUKnots   = unique(uKnot);
uniqueVKnots   = unique(vKnot);
noElemsU       = length(uniqueUKnots)-1; % # of elements xi dir.
noElemsV       = length(uniqueVKnots)-1; % # of elements eta dir.
noPtsX         = length(uKnot)-p-1;
noPtsY         = length(vKnot)-q-1;
weights        = reshape(nurbs.coefs(4,:,:),noPtsX*noPtsY,1);

controlPts = [];

for iy=1:noPtsY
  controlPts = [controlPts; nurbs.coefs(1:2,:,iy)'];
end

controlPts(:,1) = controlPts(:,1)./weights;
controlPts(:,2) = controlPts(:,2)./weights;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% chan =
%        1 2 3 4
%        5 6 7 8
%        9 10 11 12
% for a 4x3 control points
chan  = zeros(noPtsY,noPtsX);
count = 1;

for i=1:noPtsY
  for j=1:noPtsX
    chan(i,j) = count;
    count = count + 1;
  end
end

% determine our element ranges and the corresponding
% knot indices along each direction

[~,elConnU] = buildConnectivity(p,uKnot,noElemsU);
[~,elConnV] = buildConnectivity(q,vKnot,noElemsV);

% combine info from two directions to build the elements
% element is numbered as follows
%  5 | 6 | 7 | 8
% ---------------
%  1 | 2 | 3 | 4
% for a 4x2 mesh

noElems = noElemsU * noElemsV;
element = zeros(noElems,(p+1)*(q+1));

e = 1;
for v=1:noElemsV
  vConn = elConnV(v,:);
  for u=1:noElemsU
    c = 1;
    uConn = elConnU(u,:);
    for i=1:length(vConn)
      for j=1:length(uConn)
        element(e,c) = chan(vConn(i),uConn(j));
        c = c + 1;
      end
    end
    e = e + 1;
  end
end

%
[C,Cxi,Cet]  = bezierExtraction2D(uKnot,vKnot,p,q);
mesh.element = element;
mesh.node    = controlPts;
mesh.bezierExtractor = C;
mesh.p         = p;
mesh.q         = q;
mesh.elemCount = noElems;
mesh.nodeCount = noPtsX*noPtsY;
mesh.weights   = weights;
mesh.noPtsX    = noPtsX;
mesh.noPtsY    = noPtsY;




