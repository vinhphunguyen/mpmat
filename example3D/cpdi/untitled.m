a=0;
x=-5:0.001:5;
l0=[0.5 1.0 1.5 2.0];
pointCount=length(x);
vals=zeros(pointCount,length(l0));
for i=1:length(l0)
    vals(:,i) = exp(-abs(x-a)/l0(i));
end

figure
set(gca,'FontSize',14)
hold on
plot(x,vals(:,1),'black-','LineWidth',2);
plot(x,vals(:,2),'blue-','LineWidth',2);
plot(x,vals(:,3),'red-','LineWidth',2);
plot(x,vals(:,4),'cyan-','LineWidth',2);
xlabel('x')
ylabel('\phi(x)')
legend('l_0=0.5','l_0=1.0','l_0=1.5','l_0=2.0')
set(gca,'XTick',[-5 -4 -3 -2. -1 0 1 2 3 4 5])
grid on
axis([-5 5. 0 1.])
