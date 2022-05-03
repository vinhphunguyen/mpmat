%
% One dimensional Moving Least Square approximation examples.
% Convergence test.
%
% Vinh Phu Nguyen
% The University of Adelaide, Australia
% 9 October 2015.

addpath ../grid/

% node generation
L       = 1;
nn      = [20 40 80 160 320 640 1280];
ne      = nn(7);
[mesh]=buildGrid1D(L,ne,0);
node    = mesh.node;
numnode = length(node);

% MLS weight functions
shape = 'circle' ;         % shape of domain of influence
dmax  = 3.0 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;   % using cubic spline weight function

% function to be constructed using MLS, sin(pi*x)
xx     = 0:0.01:L;
fx     = sin(pi*xx);

% sampling points (think of material points in MPM or data points in
% experiments)

numMP = 2;
xp    = [];
for e=1:mesh.elemCount
  sctr  = mesh.element(e,:);
  xnode = node(sctr);
  xxp   = xnode(1) + (xnode(2)-xnode(1)).*rand(numMP,1);
  xp    = [xp;xxp];
end
up    = sin(pi*xp);

% Domain of influence for every nodes(sampling points) 
% Uniformly distributed nodes
% Definition : rad = dmax*deltaX
delta  = mesh.deltax;
di     = ones(length(xp),1)*dmax*delta;

%% Now, the real MLS stuff
% Note that data sampling points: xp,up
ui    = zeros(numnode,1);
for i=1:numnode
  pt    = node(i);
  index = defineSupport(xp,pt,di);
  if length(index) <= 2, disp('A singular'); end
  phi   = mlsLinearBasis1D(pt,index,xp,di,form);
  ui(i) = dot(phi,up(index));
end

%% compute error L2 

error = 0;
for i=1:numnode
  pt    = node(i);
  uex   = sin(pi*pt);
  error = error + (ui(i) - uex)^2;
end

error = sqrt(error/numnode);

%% visualisation

h  = [0.05;0.025;0.0125;0.00625;.003125;.0015625;.00078125];
l2 = [0.006496348257082;
      0.001680080521108;
      4.205505228667329e-04;
      1.046874641988348e-04;
      2.6232e-05;
      6.5581e-06;
      1.6419e-06];

polyfit(log(h),log(l2),1)

loglog(h,l2,'reds-','LineWidth',1.8)
xlabel('Element size h')
ylabel('Error')
set(gca,'FontSize',16)
grid on

%%
figure
hold on
plot(xx,fx,'blue-','LineWidth',1.8);
plot(node,ui,'blacks','MarkerSize',9);
plot(xp,0,'rs','MarkerSize',9,'MarkerFaceColor','blue');
xlabel('x')
ylabel('u(x)')
set(gca,'FontSize',16)
set(gca,'YTick',[0.0 0.2 0.4 0.6 0.8 1.0])
grid on
legend('u=sin(\pix)', 'MLS node values','sampling points')
%set(gca, 'YTickLabel', num2str(get(gca, 'YTick')))

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

%% visualise different weight functions

alpha = 0.3;
ra       = -1:0.01:1;
wcubic   = zeros(length(ra),1);
wquartic = zeros(length(ra),1);
wexpo    = zeros(length(ra),1);

for i=1:length(ra);
 r     = abs(ra(i));
 [wcubic(i),~]   = cubicSpline(r);
 [wquartic(i),~] = quarticSpline(r);
 [wexpo(i),~]    = expSpline(r,alpha);
end

wcubic(wcubic<0)=0;
wquartic(wquartic<0)=0;
wexpo(wexpo<0)=0;

figure
hold on
plot(ra,wcubic,'black-','LineWidth',1.8);
plot(ra,wquartic,'blue-','LineWidth',1.8);
plot(ra,wexpo,'red-','LineWidth',1.8);
xlabel('r')
ylabel('w(r)')
set(gca,'FontSize',16)
set(gca,'YTick',[0.0 0.2 0.4 0.6 0.8 1.0])
grid on
legend('cubic', 'quartic','exponential,\alpha=0.3')

%%

wexpo0    = zeros(length(ra),1);
wexpo1    = zeros(length(ra),1);
wexpo2    = zeros(length(ra),1);
wexpo3    = zeros(length(ra),1);

for i=1:length(ra);
 r     = abs(ra(i));
 [wexpo0(i),~]    = expSpline(r,0.1);
 [wexpo1(i),~]    = expSpline(r,0.2);
 [wexpo2(i),~]    = expSpline(r,0.6);
 [wexpo3(i),~]    = expSpline(r,0.8);
end

wexpo(wexpo1<0)=0;
wexpo(wexpo2<0)=0;
wexpo(wexpo3<0)=0;

figure
hold on
plot(ra,wexpo0,'cyan-','LineWidth',1.8);
plot(ra,wexpo1,'red-','LineWidth',1.8);
plot(ra,wexpo2,'blue-','LineWidth',1.8);
plot(ra,wexpo3,'black-','LineWidth',1.8);
xlabel('r')
ylabel('w(r)')
set(gca,'FontSize',16)
set(gca,'YTick',[0.0 0.2 0.4 0.6 0.8 1.0])
grid on
legend('\alpha=0,2', '\alpha=0.6','\alpha=0.8')
