A      = imread('rve.png');
grayIm = rgb2gray(A);

[noPixelX,noPixelY] = size(grayIm);

dx     = 10/noPixelX;
dy     = 10/noPixelY;

pts1 = [];
pts2 = [];

for j=1:noPixelY
    for i=1:noPixelX
       intensity = grayIm(i,j); 
       x  = (i-1)*dx + dx/2; 
       y  = (j-1)*dy + dy/2; 
       
       if     (intensity==0)
           pts1 = [pts1;x y];
       elseif (intensity==45944)
           pts2 = [pts2;x y];
       end
    end
end

figure(1)
hold on
plot(pts1(:,1),pts1(:,2),'rs');
plot(pts2(:,1),pts2(:,2),'b.');
axis equal

