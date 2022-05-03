%% Plot FEM/MPM/GIMP/CPDI basis functions
% 1D and 2D
%
% VP Nguyen
% Cardiff University, February 2014

%%
addpath ../grid/
addpath ../basis/
addpath ../particleGen/
addpath ../constitutiveModels/
addpath ../util/
addpath ../geoMesh/
addpath ../externals/
addpath ../postProcessing/

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);

%% FEM 1D function
L         = 4;
elemCount = 4;
nodes     = linspace(0,L,elemCount+1);
elements  = zeros(elemCount,2);
%
ptsCount   = 20;
xi         = -1:1/ptsCount:1;
pos        = zeros(elemCount  ,length(xi));
shapeFunc  = zeros(elemCount*2,length(xi));
dshapeFunc = zeros(elemCount*2,length(xi));

for e=1:elemCount
    pts = nodes([e e+1]);
    for i=1:length(xi) 
      pt = xi(i);
      [N,dNdxi] = lagrange_basis('L2',pt); 
      shapeFunc(2*e-1,i) = N(1);
      shapeFunc(2*e  ,i) = N(2);
      pos(e,i)           = dot(N,pts);
    end
end

%%
figure(1)
hold on
plot(pos(1,:),shapeFunc(1,:),'black-','LineWidth',1.4);
plot(pos(1,:),shapeFunc(2,:),'blue-','LineWidth',1.4);
plot(pos(2,:),shapeFunc(3,:),'blue-','LineWidth',1.4);
plot(pos(2,:),shapeFunc(4,:),'red-','LineWidth',1.4);
plot(pos(3,:),shapeFunc(5,:),'red-','LineWidth',1.4);
plot(pos(3,:),shapeFunc(6,:),'green-','LineWidth',1.4);
plot(pos(4,:),shapeFunc(7,:),'green-','LineWidth',1.4);
plot(pos(4,:),shapeFunc(8,:),'cyan-','LineWidth',1.4);
%plot(x,sum(dshapeFunc,1),'black--' ,'LineWidth',1.4);
axis equal
axis([0 4 0 1])
set(gca,'FontSize',14)
grid on, box on
% exportfig(gcf,'L2-basis-example.eps',opts)

%% FEM 1D function (quadratic)
L         = 4;
elemCount = 4;
nodes     = linspace(0,L,2*elemCount+1);
elements  = zeros(elemCount,2);
%
ptsCount   = 20;
xi         = -1:1/ptsCount:1;
pos        = zeros(elemCount  ,length(xi));
shapeFunc  = zeros(elemCount*3,length(xi));
dshapeFunc = zeros(elemCount*3,length(xi));

for e=1:elemCount
    pts = nodes([2*e-1 2*e 2*e+1]);
    for i=1:length(xi) 
      pt = xi(i);
      N=[(1-pt)*pt/(-2);1-pt^2;(1+pt)*pt/2];
      shapeFunc(3*e-2,i) = N(1);
      shapeFunc(3*e-1,i) = N(3);
      shapeFunc(3*e  ,i) = N(2);
      pos(e,i)           = dot(N,pts);
    end
end

%%
figure(2)
hold on
plot(pos(1,:),shapeFunc(1,:),'black-','LineWidth',1.4);
plot(pos(1,:),shapeFunc(2,:),'blue-','LineWidth',1.4);
plot(pos(1,:),shapeFunc(3,:),'red-','LineWidth',1.4);
plot(pos(2,:),shapeFunc(4,:),'black-','LineWidth',1.4);
plot(pos(2,:),shapeFunc(5,:),'blue-','LineWidth',1.4);
plot(pos(2,:),shapeFunc(6,:),'red-','LineWidth',1.4);
plot(pos(3,:),shapeFunc(7,:),'black-','LineWidth',1.4);
plot(pos(3,:),shapeFunc(8,:),'blue-','LineWidth',1.4);
plot(pos(3,:),shapeFunc(9,:),'red-','LineWidth',1.4);
plot(pos(4,:),shapeFunc(10,:),'black-','LineWidth',1.4);
plot(pos(4,:),shapeFunc(11,:),'blue-','LineWidth',1.4);
plot(pos(4,:),shapeFunc(12,:),'red-','LineWidth',1.4);
%plot(x,sum(dshapeFunc,1),'black--' ,'LineWidth',1.4);
axis equal
axis([0 4 -0.5 1])
set(gca,'FontSize',14)
grid on, box on
% exportfig(gcf,'L3-basis-example.eps',opts)

%% GIMP function
% L         = 4;
% elemCount = 4;
% nodes     = linspace(0,L,elemCount+1);
% elements  = zeros(elemCount,2);
% 
% h         = L/elemCount;
% n         = 2;
% lp        = h/n;  % assume "n" particles per cell
% 
% lpa = [3 1];
% 
% %
% 
% ptsCount   = 250;
% x          = 0:L/ptsCount:L;
% shapeFunc  = zeros(elemCount+1,length(x));
% dshapeFunc = zeros(elemCount+1,length(x));
% 
% for I=1:elemCount+1
%     for i=1:length(x)           
%         pts = x(i) - nodes(I);
%         [shapeFunc(I,i), dshapeFunc(I,i)] = getGIMP(pts,h,1.5);
%     end
% end
% 
% %%
% figure(1)
% hold on
% plot(x,shapeFunc(1,:),'black-','LineWidth',1.4);
% plot(x,shapeFunc(2,:),'red-'  ,'LineWidth',1.4);
% plot(x,shapeFunc(3,:),'blue-' ,'LineWidth',1.4);
% plot(x,shapeFunc(4,:),'cyan-' ,'LineWidth',1.4);
% plot(x,shapeFunc(5,:),'yellow-' ,'LineWidth',1.4);
% %plot(x,sum(dshapeFunc,1),'black--' ,'LineWidth',1.4);
% axis equal
% axis([0 4 0 1])
% %exportfig(gcf,'gimp-basis-example1.eps',opts)

%% MPM function

% for I=1:5
%     for i=1:length(x)           
%         pts = x(i) - nodes(I);
%         [shapeFunc(I,i), dshapeFunc(I,i)] = getMPM(pts,h);
%     end
% end
% 
% %%
% figure(10)
% hold on
% plot(x,shapeFunc(1,:),'black-','LineWidth',1.4);
% plot(x,shapeFunc(2,:),'red-'  ,'LineWidth',1.4);
% plot(x,shapeFunc(3,:),'blue-' ,'LineWidth',1.4);
% plot(x,shapeFunc(4,:),'cyan-' ,'LineWidth',1.4);
% plot(x,shapeFunc(5,:),'yellow-' ,'LineWidth',1.4);
% axis equal
% axis([0 4 0 1])

%%
% figure(2)
% hold on
% 
% plot(x,dshapeFunc(3,:),'red-'  ,'LineWidth',1.4);
% 
% 
% axis equal
% %exportfig(gcf,'gimp-basis-example2.eps',opts)
% 
% %% plot 2D GIMP functions
% 
% [X,Y] = meshgrid (linspace(0,4,200));
% R     = zeros(size(X,1),size(X,1));
% 
% xI    = [2 2];
% lp    = 1;
% h     = 1;
% 
% for i=1:size(X,1)
%     for j=1:size(X,1)
%         x  = X(1,i);
%         y = Y(j,1);        
%         pts = x - xI(1);
%         [Nx, dshape] = getGIMP(pts,h,lp);
%         pts = y - xI(2);
%         [Ny, dshape] = getGIMP(pts,h,lp);
%         R(i,j) = Nx*Ny;
%     end
% end
% 
% %exportfig(gcf,'cubic-N3',opts)
% 
% figure
% surf (X,Y,R,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
% hold off
% 
% figure
% surf (X,Y,R,'EdgeColor','none','LineStyle','none')
% hold off
% axis equal
% view([0 90])
% colorbar
% 
% 
% 
% %% plot 2D MPM functions
% 
% [X,Y] = meshgrid (linspace(0,4,200));
% R     = zeros(size(X,1),size(X,1));
% 
% xI    = [2 2];
% lp    = 1;
% h     = 1;
% 
% for i=1:size(X,1)
%     for j=1:size(X,1)
%         x  = X(1,i);
%         y = Y(j,1);        
%         pts = x - xI(1);
%         [Nx, dshape] = getMPM(pts,h);
%         pts = y - xI(2);
%         [Ny, dshape] = getMPM(pts,h);
%         R(i,j) = Nx*Ny;
%     end
% end
% 
% %exportfig(gcf,'cubic-N3',opts)
% 
% figure
% surf (X,Y,R,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
% hold off
% 
% figure
% surf (X,Y,R,'EdgeColor','none','LineStyle','none')
% hold off
% axis equal
% view([0 90])
% colorbar

% %% plot quadratic Lagrange basis functions
%
% % x=-1:0.01:1;
% % N1=0.5*x.*(x-1);
% % N2=1-x.^2;
% % N3=0.5*x.*(x+1);
% %
% % figure(3)
% % hold on
% % plot(x,N1,'black-','LineWidth',1.4);
% % plot(x,N2,'blue-','LineWidth',1.4);
% % plot(x,N3,'red-','LineWidth',1.4);
%
%
% for I=1:elemCount
%     for i=1:length(x)
%         pts = x(i) - nodes(I);
%         r = abs(pts/h);
%         dd = sign(pts)/h;
%         if r <= 0.25
%             shapeFunc(I,i) =  (7-16*r^2)/8;
%             dshapeFunc(I,i)= dd*(-4*r);
%         elseif r <= 0.75
%             shapeFunc(I,i) =  1-r;
%             dshapeFunc(I,i)= -dd;
%         elseif r <=1.25
%             shapeFunc(I,i) =  (5-4*r)^2/16;
%             dshapeFunc(I,i)= -dd*(5-4*r)/2;
%         else
%             shapeFunc(I,i) =  0;
%             dshapeFunc(I,i) =  0;
%         end
%     end
% end
%
% figure(10)
% hold on
% plot(x,shapeFunc(1,:),'black-','LineWidth',1.4);
% plot(x,shapeFunc(2,:),'red-'  ,'LineWidth',1.4);
% plot(x,shapeFunc(3,:),'blue-' ,'LineWidth',1.4);
% plot(x,shapeFunc(4,:),'cyan-' ,'LineWidth',1.4);
% plot(x,shapeFunc(5,:),'yellow-' ,'LineWidth',1.4);
% plot(x,sum(shapeFunc,1),'black--' ,'LineWidth',1.4);
% plot([1 1],[-0.1 1.1],'black-' ,'LineWidth',1.8);
% plot([3 3],[-0.1 1.1],'black-' ,'LineWidth',1.8);
% axis equal
% %axis([0 4 -0.1 1.1])
%
% figure(20)
% hold on
% plot(x,dshapeFunc(1,:),'black-','LineWidth',1.4);
% plot(x,dshapeFunc(2,:),'red-'  ,'LineWidth',1.4);
% plot(x,dshapeFunc(3,:),'blue-' ,'LineWidth',1.4);
% plot(x,dshapeFunc(4,:),'cyan-' ,'LineWidth',1.4);
% plot(x,dshapeFunc(5,:),'yellow-' ,'LineWidth',1.4);
% %plot(x,sum(dshapeFunc,1),'black--' ,'LineWidth',1.4);
% axis([0 4 -1 1])
% axis equal
%
%
% shapeFunc  = zeros(elemCount,length(x));
% dshapeFunc = zeros(elemCount,length(x));
%
% for I=1:elemCount
%     xI = 0.5*(nodes(I)+nodes(I+1));
%     for i=1:length(x)
%         pts = x(i) - xI;
%         [shapeFunc(I,i), dshapeFunc(I,i)] = getQuadraticBspline(pts,h);
%     end
% end
%
%
% figure(20)
% hold on
% plot(x,shapeFunc(1,:),'black-','LineWidth',1.4);
% plot(x,shapeFunc(2,:),'red-'  ,'LineWidth',1.4);
% plot(x,shapeFunc(3,:),'blue-' ,'LineWidth',1.4);
% plot(x,shapeFunc(4,:),'cyan-' ,'LineWidth',1.4);
% plot(x,shapeFunc(5,:),'yellow-' ,'LineWidth',1.4);
% %plot(x,sum(shapeFunc,1),'black--' ,'LineWidth',1.4);
% %axis equal

%%
%%%%%%%%%%%%
%% plot CPDI functions 1D (linear)

% grid
% L          = 3;
% elemCount  = 3;
% dx         = L/elemCount;
% nodes      = linspace(0,L,elemCount+1);
% % particle
% lp         = 1.5;
% % phi_I = 0.5*N_I(x1) +0.5*N_I(x2)
% ptsCount   = 250;
% xx          = 0:L/ptsCount:L;
% x          = 0:L/ptsCount:L;
% shapeFunc  = zeros(elemCount+1,length(x));
% shapeFuncS = zeros(elemCount+1,length(x));
% dshapeFunc = zeros(elemCount+1,length(x));
% 
% for I=1:elemCount+1
%     for i=1:length(x)        
%         xc1 = (x(i) - 0.5*lp) - nodes(I); N1 = getMPM(xc1,dx);
%         xc2 = (x(i) + 0.5*lp) - nodes(I); N2 = getMPM(xc2,dx);
%         shapeFunc(I,i)  = 0.5*N1 + 0.5*N2;        
%         dshapeFunc(I,i) = (-1/lp)*N1 + (1/lp)*N2;    
%     end
% end
% 
% for I=1:elemCount+1
%     for i=1:length(xx)        
%         shapeFuncS(I,i) = getMPM(xx(i)-nodes(I),h);
%     end
% end
% 
% figure(2)
% 
% subplot(2,1,1)
% set(gca,'FontSize',14)
% hold on
% plot(xx,shapeFuncS(1,:),'black-','LineWidth',1.8);
% plot(xx,shapeFuncS(2,:),'red-'  ,'LineWidth',1.8);
% plot(xx,shapeFuncS(3,:),'blue-' ,'LineWidth',1.8);
% plot(xx,shapeFuncS(4,:),'cyan-' ,'LineWidth',1.8);
% 
% subplot(2,1,2)
% set(gca,'FontSize',14)
% hold on
% plot(x,shapeFunc(1,:),'black-','LineWidth',1.8);
% plot(x,shapeFunc(2,:),'red-'  ,'LineWidth',1.8);
% plot(x,shapeFunc(3,:),'blue-' ,'LineWidth',1.8);
% plot(x,shapeFunc(4,:),'cyan-' ,'LineWidth',1.8);
% plot(x,shapeFunc(1,:)+shapeFunc(2,:)+shapeFunc(3,:)+shapeFunc(4,:),'m-' ,'LineWidth',1.8);

%exportfig(gcf,'splinecurve.eps',opts)
% 
% figure(3)
% 
% subplot(2,1,1)
% set(gca,'FontSize',14)
% hold on
% plot(xx,shapeFunc(1,:),'black-','LineWidth',1.8);
% plot(xx,shapeFunc(2,:),'red-'  ,'LineWidth',1.8);
% plot(xx,shapeFunc(3,:),'blue-' ,'LineWidth',1.8);
% plot(xx,shapeFunc(4,:),'cyan-' ,'LineWidth',1.8);
% 
% subplot(2,1,2)
% set(gca,'FontSize',14)
% hold on
% plot(x,dshapeFunc(1,:),'black-','LineWidth',1.8);
% plot(x,dshapeFunc(2,:),'red-'  ,'LineWidth',1.8);
% plot(x,dshapeFunc(3,:),'blue-' ,'LineWidth',1.8);
% plot(x,dshapeFunc(4,:),'cyan-' ,'LineWidth',1.8);
%plot(x,dshapeFunc(1,:)+dshapeFunc(2,:)+dshapeFunc(3,:)+dshapeFunc(4,:),'m-' ,'LineWidth',1.8);

%plot(x,sum(dshapeFunc,1),'black--' ,'LineWidth',1.4);
%axis equal
%axis([0 4 0 1])

%         pts = x(i) - nodes(I);
%         [shapeFunc(I,i), dshapeFunc(I,i)] = getMPM(pts,h);

%%%%%%%%%%%%
%% plot CPDI functions 1D (quadratic)

% grid
% L          = 4;
% elemCount  = 4;
% dx         = L/elemCount;
% nodes      = linspace(0,L,elemCount+1);
% % particle
% lp         = 1;
% % phi_I = 1/6*(x21-2c)*N_I(x1) + 1/6*(x_21+2c)*N_I(x2) + x21 N_I(x3)
% % x21=x2-x1 and  c = x1+x2-2x3
% ptsCount   = 250;
% xx         = 0:L/ptsCount:L;
% x          = 0:L/ptsCount:L;
% shapeFunc  = zeros(elemCount+1,length(x));
% shapeFuncS = zeros(elemCount+1,length(x));
% dshapeFunc = zeros(elemCount+1,length(x));
% 
% for I=1:elemCount+1
%     for i=1:length(x)        
%         x1 = x(i) - 0.5*lp; 
%         x2 = x(i) + 0.5*lp; 
%         x3 = x(i) ; 
%         x21= x2-x1;
%         c  = x1+x2-2*x3;
%         wf = [(1/6)*(x21-2*c) (1/6)*(x21+2*c) 4*x21/6];
%         wg = [-1 1 0];
%         Na = [getMPM(x1- nodes(I),dx);getMPM(x2- nodes(I),dx);getMPM(x3- nodes(I),dx)];
%         shapeFunc(I,i)  = (1/lp)*dot(wf,Na);        
%         dshapeFunc(I,i) = (1/lp)*dot(wg,Na);        
%     end
% end
% 
% for I=1:elemCount+1
%     for i=1:length(xx)        
%         shapeFuncS(I,i) = getMPM(xx(i)-nodes(I),dx);
%     end
% end
% 
% figure(2)
% 
% subplot(2,1,1)
% set(gca,'FontSize',14)
% hold on
% for I=1:elemCount+1
% plot(xx,shapeFuncS(I,:),'LineWidth',2.1);
% end
% 
% subplot(2,1,2)
% set(gca,'FontSize',14)
% hold on
% for I=1:elemCount+1
% plot(xx,shapeFunc(I,:),'LineWidth',2.1);
% end
% plot(x,sum(shapeFunc(:,:)),'m-' ,'LineWidth',2.1);
% 
% exportfig(gcf,'splinecurve.eps',opts)
% 
% figure(3)
% 
% subplot(2,1,1)
% set(gca,'FontSize',14)
% hold on
% for I=1:elemCount+1
% plot(xx,shapeFunc(I,:),'LineWidth',2.1);
% end
% 
% subplot(2,1,2)
% set(gca,'FontSize',14)
% hold on
% for I=1:elemCount+1
% plot(x,dshapeFunc(I,:),'LineWidth',2.1);
% end
%%

%% plot CPDI functions 1D (linear with subsampling)

% grid
% L          = 4;
% elemCount  = 4;
% dx         = L/elemCount;
% nodes      = linspace(0,L,elemCount+1);
% % particle
% lp         = 1;
% % phi_I = 1/6*(x21-2c)*N_I(x1) + 1/6*(x_21+2c)*N_I(x2) + x21 N_I(x3)
% % x21=x2-x1 and  c = x1+x2-2x3
% ptsCount   = 250;
% xx         = 0:L/ptsCount:L;
% x          = 0:L/ptsCount:L;
% shapeFunc  = zeros(elemCount+1,length(x));
% shapeFuncS = zeros(elemCount+1,length(x));
% dshapeFunc = zeros(elemCount+1,length(x));
% 
% for I=1:elemCount+1
%     for i=1:length(x)        
%         x1 = x(i) - 0.5*lp; 
%         x2 = x(i) + 0.5*lp; 
%         x3 = x(i) ; 
%         wf = [1/2 1/2 1/2];
%         wg = [-1 1 0];
%         Na = [getMPM(x1- nodes(I),dx);getMPM(x2- nodes(I),dx);getMPM(x3- nodes(I),dx)];
%         shapeFunc(I,i)  = (1/lp)*dot(wf,Na);        
%         dshapeFunc(I,i) = (1/lp)*dot(wg,Na);        
%     end
% end
% 
% for I=1:elemCount+1
%     for i=1:length(xx)        
%         shapeFuncS(I,i) = getMPM(xx(i)-nodes(I),dx);
%     end
% end
% 
% figure(2)
% 
% subplot(2,1,1)
% set(gca,'FontSize',14)
% hold on
% for I=1:elemCount+1
% plot(xx,shapeFuncS(I,:),'LineWidth',2.1);
% end
% 
% subplot(2,1,2)
% set(gca,'FontSize',14)
% hold on
% for I=1:elemCount+1
% plot(xx,shapeFunc(I,:),'LineWidth',2.1);
% end
% plot(x,sum(shapeFunc(:,:)),'m-' ,'LineWidth',2.1);
% 
% exportfig(gcf,'splinecurve.eps',opts)
% 
% figure(3)
% 
% subplot(2,1,1)
% set(gca,'FontSize',14)
% hold on
% for I=1:elemCount+1
% plot(xx,shapeFunc(I,:),'LineWidth',2.1);
% end
% 
% subplot(2,1,2)
% set(gca,'FontSize',14)
% hold on
% for I=1:elemCount+1
% plot(x,dshapeFunc(I,:),'LineWidth',2.1);
% end
