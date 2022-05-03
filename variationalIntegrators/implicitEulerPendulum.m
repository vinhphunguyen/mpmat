% Implicit EUler integration for the motion of the pendulum problem

addpath ../util/

% constants

g = 1; l = 1; m = 1;

t0 = 0;
tN = 75; % seconds

h  = 0.01;
N  = tN/h;

time = t0:h:tN;

q  = zeros(N+1,1);
v  = zeros(N+1,1);
e  = zeros(N+1,1); % total energy

% initial conditions

q(1) = 1.0; 
v(1) = 1.0; 
e(1) = 0.5*v(1)^2 - m*g*l*cos(q(1));

% time integration
a      = h*h*g/l;
eps    = 1e-7;


for k=1:N
    qk = q(k);
    vk = v(k);
    % solve for q_{k+1} = x
    x0 = 0; deltaq = 1.;
    while ( abs(deltaq) > eps )
        f      = x0 + a*sin(x0) - qk - h*vk;
        deltaq = -f/(1+a*cos(x0));
        x0     = x0 + deltaq;
    end
    q(k+1) = x0;
    v(k+1) = vk - (h*g/l)*sin(x0);
    
    e(k+1) = 0.5*v(k+1)^2 - m*g*l*cos(q(k+1));
end

%%
opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);


figure(1)
set(gca,'FontSize',14)
hold on
plot(time,q,'b-','LineWidth',2.2);
xlabel('Time [s]')
ylabel('q(t)')
grid on

%%
figure(2)
set(gca,'FontSize',14)
hold on
plot(time,v,'b-','LineWidth',2.2);
xlabel('Time [s]')
ylabel('q(t)')
grid on

%%

figure(3)
set(gca,'FontSize',14)
hold on
plot(q(1),v(1),'o','MarkerSize',20,'MarkerFaceColor','red')
plot(q,v,'b-','LineWidth',2.2);
xlabel('q(t)')
ylabel('v(t)')
grid on

%%

figure(4)
set(gca,'FontSize',14)
hold on
plot(time,e,'b-','LineWidth',2.2);
xlabel('time [s]')
ylabel('energy')
grid on