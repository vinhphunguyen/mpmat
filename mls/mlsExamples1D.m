%
% One dimensional Moving Least Square approximation examples.
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 9 October 2015.

% node generation
L       = 1;
nnx     = 11;
node    = linspace(0,L,nnx);
numnode = length(node);

% MLS weight functions
shape = 'circle' ;         % shape of domain of influence
dmax  = 2.5 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;   % using cubic spline weight function

% function to be constructed using MLS
xx     = 0:0.01:L;
fx     = 1 + xx;                     % linear function
fxc    = ones(length(xx),1);         % constant function

% sampling points (think of material points in MPM or data points in
% experiments)

numMP = 15;
xp    = 0 + (1-0).*rand(numMP,1);
up    = 1 + xp;
upc   = ones(length(xp),1);

% Domain of influence for every nodes(sampling points) 
% Uniformly distributed nodes
% Definition : rad = dmax*deltaX
delta  = L/(nnx-1);
di     = ones(numMP,1)*dmax*delta;

%% Now, the real MLS stuff
% Note that data sampling points: xp,up
ui1    = zeros(numnode,1);
ui2    = zeros(numnode,1);
for i=1:numnode
  pt    = node(i);
  index = defineSupport(xp,pt,di);
  phi   = mlsLinearBasis1D(pt,index,xp,di,form);
  ui1(i) = dot(phi,up(index));
  phi   = mlsConstantBasis1D(pt,index,xp,di,form);  
  ui2(i) = dot(phi,up(index));
end


%% visualisation
figure(1)
hold on
plot(xx,fx,'blue-','LineWidth',1.8);
plot(node,ui1,'blacks','MarkerSize',9);
plot(xp,1,'rs','MarkerSize',9,'MarkerFaceColor','blue');
xlabel('x')
ylabel('u(x)')
set(gca,'FontSize',16)
set(gca,'YTick',[1.0 1.2 1.4 1.6 1.8 2.0])
grid on
legend('u=1+x', 'MLS node values','sampling points')
%set(gca, 'YTickLabel', num2str(get(gca, 'YTick')))
title('Linear MLS')

figure(2)
hold on
plot(xx,fx,'blue-','LineWidth',1.8);
plot(node,ui2,'blacks','MarkerSize',9);
plot(xp,1,'rs','MarkerSize',9,'MarkerFaceColor','blue');
xlabel('x')
ylabel('u(x)')
set(gca,'FontSize',16)
set(gca,'YTick',[1.0 1.2 1.4 1.6 1.8 2.0])
grid on
legend('u=1+x', 'MLS node values','sampling points')
%set(gca, 'YTickLabel', num2str(get(gca, 'YTick')))
title('Zero-order MLS')

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

%% Now, the real MLS stuff, constant function
% Note that data sampling points: xp,upc
ui1    = zeros(numnode,1);
for i=1:numnode
  pt    = node(i);
  index = defineSupport(xp,pt,di);
  phi   = mlsLinearBasis1D(pt,index,xp,di,form);
  ui1(i) = dot(phi,upc(index));
end


%% visualisation
figure(3)
hold on
plot(xx,fxc,'blue-','LineWidth',1.8);
plot(node,ui1,'blacks','MarkerSize',9);
plot(xp,1,'rs','MarkerSize',9,'MarkerFaceColor','blue');
xlabel('x')
ylabel('u(x)')
set(gca,'FontSize',16)
set(gca,'YTick',[1.0 1.2 1.4 1.6 1.8 2.0])
grid on
legend('u=1', 'MLS node values','sampling points')
%set(gca, 'YTickLabel', num2str(get(gca, 'YTick')))
