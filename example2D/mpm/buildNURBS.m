function nurbs = buildNURBS(knots,controlPts,weights,noPtsX,noPtsY)


controlPts(:,1) = controlPts(:,1).*weights;
controlPts(:,2) = controlPts(:,2).*weights;

pts        = zeros(4,noPtsX,noPtsY);
pts(1,:,:) = reshape(controlPts(:,1),noPtsX,noPtsY);
pts(2,:,:) = reshape(controlPts(:,2),noPtsX,noPtsY);
pts(4,:,:) = reshape(weights,noPtsX,noPtsY);

nurbs = nrbmak(pts,knots);
