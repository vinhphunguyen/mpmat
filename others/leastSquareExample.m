n  = 6;
x  = [0;0.5;1;1.5;2;2.5];
u  = x.^2;

m  = 3;

A  = zeros(m,m);
b  = zeros(m,1);

for i=1:n
  xi = x(i);
  p  = [1;xi;xi^2];
  A  = A + p*p';
  b  = b+ p*u(i);
end

a    = A\b;

%%

xx = 0:0.01:x(end);
uh = a(1) + a(2)*xx + a(3)*xx.^2;

figure
hold on
plot(x,u,'rs','MarkerSize',9,'MarkerFaceColor','red')
plot(xx,uh,'b-','LineWidth',2)
set(gca,'FontSize',16)
legend('data','least-square fit')
grid on

%% sin function

x  = -2:0.25:4;
u  = sin(x);

n  = length(x);
m  = 3;

A  = zeros(m,m);
b  = zeros(m,1);

for i=1:n
  xi = x(i);
  p  = [1;xi;xi^2];
  A  = A + p*p';
  b  = b+ p*u(i);
end

a1    = A\b;

m  = 4;

A  = zeros(m,m);
b  = zeros(m,1);

for i=1:n
  xi = x(i);
  p  = [1;xi;xi^2;xi^3];
  A  = A + p*p';
  b  = b+ p*u(i);
end

a2    = A\b;

m  = 5;

A  = zeros(m,m);
b  = zeros(m,1);

for i=1:n
  xi = x(i);
  p  = [1;xi;xi^2;xi^3;xi^4];
  A  = A + p*p';
  b  = b+ p*u(i);
end

a3    = A\b;

%%

xx = -2:0.01:4;
uh2 = a1(1) + a1(2)*xx + a1(3)*xx.^2;
uh3 = a2(1) + a2(2)*xx + a2(3)*xx.^2 + a2(4)*xx.^3;
uh4 = a3(1) + a3(2)*xx + a3(3)*xx.^2 + a3(4)*xx.^3 + a3(5)*xx.^4;

figure
hold on
plot(x,u,'rs','MarkerSize',9,'MarkerFaceColor','red')
plot(xx,uh2,'b-','LineWidth',2)
plot(xx,uh3,'cyan-','LineWidth',2)
plot(xx,uh4,'green-','LineWidth',2)
set(gca,'FontSize',16)
legend('data','least-square 2nd', 'least-square 3rd', 'least-square 4th')
grid on
