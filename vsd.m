clear; close all; clc;
%%
kp = 8000;
kv = 1.0*(2*kp)^.5;

%%
T = 0.8;
dt = 0.001;
t = (0:dt:T)';
demo_x = -(t - sin(2*pi*t.^3).*cos(2*pi*t.^3).*exp(t.^4));

xdd = zeros(length(t),1);
xd = zeros(length(t),1);
x = zeros(length(t),1);
xdd(1,1) = 0;
xd(1,1) = 0;
x(1,1) = 0;

%%
for i = 1:length(t)-1
    xdd(i+1) = kp*(demo_x(i)-x(i)) - kv*xd(i);
    xd(i+1) = xd(i) + xdd(i)*dt;
    x(i+1) = x(i) + xd(i)*dt;
    
    if  rem(i,1)==0
        %plot([t(i),t(i)], [x(i),demo_x(i)] ,'r');
        %hold on;
        %plot([x(i),x(i) + x_t], [y(i),y(i) + y_t]);
        %hold on;
    end
end

%%
plot(t, demo_x, 'b');
hold on;
plot(t, x, 'r');
grid on;
legend('demonstration','reproduction');
title('Reproduction of Demonstration with VSD');
xlabel('x','fontsize',18);
ylabel('y','fontsize',18);
% xlim([0,1]);
set(gca,'FontSize',12);

