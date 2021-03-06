clear; close all; clc;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp_1 = 500;
kv_1 = 2.3*(2*kp_1)^.5;
kp_2 = 20000; % 200000,100000,50000,20000,10000,5000,3000,1000,0
kv_2 = 1*(2*kp_2)^.5;
kp_3 = 0; % 5500
kv_3 = 0.8*(2*kp_3)^.5;
alpha = 0.01; 
Fe = 0;

%%
T = 0.4;
dt = 0.001;
t = (0:dt:T)';
demo_x = t;
demo_y = -(demo_x - sin(2*pi*demo_x.^3).*cos(2*pi*demo_x.^3).*exp(demo_x.^4));
x_end = demo_x(end) + 0.0;
y_end = demo_y(end) + 0.05;

xdd = zeros(length(t),1);
xd = zeros(length(t),1);
x = zeros(length(t),1);
xdd(1,1) = 0;
xd(1,1) = 0;
x(1,1) = 0.0;% 0.1,0.2;

ydd = zeros(length(t),1);
yd = zeros(length(t),1);
y = zeros(length(t),1);
ydd(1,1) = 0;
yd(1,1) = 0;
y(1,1) = 0.0;% 0.4,0.3,0.2,0.1,0.0;

R_search_start = 1;
R_search_end = length(t);%length(x)/10;
index_data = zeros(2,1);
index_data(1,1) = 1; 
index_data(2,1) = +inf; 
R_search = 100;
steps = 2000;
a = steps; %��Բ����
b = 1; %��Բ����

%%
for i = 1:steps %length(t)-1
    for j = R_search_start:R_search_end %length(x)
        s = (demo_x(j) - x(i))^2 + (demo_y(j) - y(i))^2;
        if s < index_data(2,1)
            index_data(2,1) = s;
            index_data(1,1) = j;
        end
        R_search_start = index_data(1,1) - R_search;
        if R_search_start < 1
            R_search_start = 1;
        end
        R_search_end = index_data(1,1) + R_search;
        if R_search_end > length(t)
            R_search_end = length(t);
        end
    end
    
    %index_data(1,1)
    y_r = demo_y(index_data(1,1))-y(i);
    x_r = demo_x(index_data(1,1))-x(i);
    
    y_t = -x_r/(y_r+0.00001);
    x_t = 1;
    t_mod = (y_t^2+x_t^2)^.5;
    y_t = y_t/t_mod;
    x_t = x_t/t_mod;
    
    y_a = y_end - y(i);
    x_a = x_end - x(i);
    k_e = exp(-alpha*i);
    
    %ydd(i+1) = kp_1*y_a - kv_1*yd(i) + k_ey*(kp_2*(demo_y(i)-y(i)) - kv_2*yd(i)) + Fe; % e+i
    %ydd(i+1) = kp_1*y_t - kv_1*yd(i) + k_ey*(kp_2*y_r - kv_2*yd(i)) + Fe; % t+r, ok
    ydd(i+1) = kp_1*y_a - kv_1*yd(i) + k_e*(kp_2*y_r - kv_2*yd(i) + kp_3*y_t + kv_3*yd(i)) + Fe; % e+r+t
    %ydd(i+1) = kp_1*(demo_y(i)-y(i)) - kv_1*yd(i) + Fe; % vsd
    yd(i+1) = yd(i) + ydd(i)*dt;
    y(i+1) = y(i) + yd(i)*dt;
    
    %xdd(i+1) = kp_1*x_a - kv_1*xd(i) + k_ex*(kp_2*(demo_x(i)-x(i)) - kv_2*xd(i)) + Fe; % e+i
    %xdd(i+1) = kp_1*x_t - kv_1*xd(i) + k_ex*(kp_2*x_r - kv_2*xd(i)) + Fe; % t+r, ok
    xdd(i+1) = kp_1*x_a - kv_1*xd(i) + k_e*(kp_2*x_r - kv_2*xd(i) + kp_3*x_t + kv_3*xd(i)) + Fe; % e+r+t
    %xdd(i+1) = kp_1*(demo_x(i)-x(i)) - kv_1*xd(i) + Fe; % vsd
    xd(i+1) = xd(i) + xdd(i)*dt;
    x(i+1) = x(i) + xd(i)*dt;
    
    index_data(2,1) = +inf;
    
    if  rem(i,10)==0
        plot([x(i),x(i) + x_r], [y(i),y(i) + y_r],'r');
        hold on;
        plot([x(i),x(i) + x_t/100], [y(i),y(i) + y_t/100]);
        hold on;
%         plot(x(i), k_e, 'r*');
%         hold on;
    end
end

%%
demo = plot(demo_x, demo_y, 'b', 'LineWidth',1);
hold on;
repro = plot(x, y, 'LineWidth',1);
hold on;
redpoint = plot(-0.1, 0.25, 'r');
hold on;
legend([demo, repro,redpoint],'demonstration','reproduction','y_\perp');
% title('Nearest points of the Demonstration with MVSD');
xlabel('x','fontsize',18);
ylabel('y','fontsize',18);
set(gca,'FontSize',12);
% axis equal;
axis off;
grid on;
