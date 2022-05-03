fileID = fopen('nuru-lodi.dat','r');
formatSpec = '%f %f %f %f';
sizeA = [4 Inf]; % column order

A = fscanf(fileID,formatSpec,sizeA);
A = A';
fclose(fileID);

A(:,2) = A(:,2)/1000;
A(:,4) = A(:,4)/1000;

figure
clf;
set( gcf, 'DoubleBuffer', 'on' )
set(gca,'FontSize',14)
hold on
plot(exp(:,1),exp(:,2),'bo-','LineWidth',1.6);
plot(A(:,3),A(:,4),'r*-','LineWidth',1.6);
xlabel('u [mm]')
ylabel('P [N]')
