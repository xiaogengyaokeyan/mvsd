clear; close all; clc;
%%
kp = 300;
kv = 1.0*(2*kp)^.5;

%%
T = 0.8;
dt = 0.002;
% t = (0:dt:T)';

load('Data/demot4');
load('Data/demox4.mat');
load('Data/demoy4.mat');
load('Data/demoz4.mat');
t = t_1;
demo_x = x_1;
demo_y = y_1;
demo_z = z_1;

demo_xd = diff(x_1)/dt;
demo_yd = diff(y_1)/dt;
demo_zd = diff(z_1)/dt;

xdd = zeros(length(t),1);
xd = zeros(length(t),1);
x = zeros(length(t),1);
ydd = zeros(length(t),1);
yd = zeros(length(t),1);
y = zeros(length(t),1);
zdd = zeros(length(t),1);
zd = zeros(length(t),1);
z = zeros(length(t),1);
xdd(1,1) = 0;
xd(1,1) = 0;
x(1,1) = demo_x(1,1);
ydd(1,1) = 0;
yd(1,1) = 0;
y(1,1) = demo_y(1,1);
zdd(1,1) = 0;
zd(1,1) = 0;
z(1,1) = demo_z(1,1);

%%
for i = 1:length(t)-1
%     xdd(i+1) = kp*(demo_x(i)-x(i)) + kv*(demo_xd(i)-xd(i));
    xdd(i+1) = kp*(demo_x(i)-x(i)) - kv*xd(i);
    xd(i+1) = xd(i) + xdd(i)*dt;
    x(i+1) = x(i) + xd(i)*dt;
    
%     ydd(i+1) = kp*(demo_y(i)-y(i)) + kv*(demo_yd(i)-yd(i));
    ydd(i+1) = kp*(demo_y(i)-y(i)) - kv*yd(i);
    yd(i+1) = yd(i) + ydd(i)*dt;
    y(i+1) = y(i) + yd(i)*dt;
    
%     zdd(i+1) = kp*(demo_z(i)-z(i)) + kv*(demo_zd(i)-zd(i));
    zdd(i+1) = kp*(demo_z(i)-z(i)) - kv*zd(i);
    zd(i+1) = zd(i) + zdd(i)*dt;
    z(i+1) = z(i) + zd(i)*dt;
    
    if  rem(i,1)==0
        %plot([t(i),t(i)], [x(i),demo_x(i)] ,'r');
        %hold on;
        %plot([x(i),x(i) + x_t], [y(i),y(i) + y_t]);
        %hold on;
    end
end

%%
figure;
plot(t,demo_x);
hold on;
plot(t,demo_y);
hold on;
plot(t,demo_z);
hold on;
plot(t,x);
hold on;
plot(t,y);
hold on;
plot(t,z);

figure;
plot(t,xdd,'r');
hold on;
plot(t,ydd,'g');
hold on;
plot(t,zdd,'b');
hold on;

figure;
plot3(demo_x,demo_y,demo_z, 'b');
hold on;
plot3(x, y, z, 'r');
grid on;
legend('demonstration','reproduction');
title('Reproduction of Demonstration with VSD\_3d');
xlabel('x','fontsize',18);
ylabel('y','fontsize',18);
zlabel('z','fontsize',18);
% xlim([0,1]);
set(gca,'FontSize',12);

