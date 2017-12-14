clear; close all; clc;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp_1 = 500;
kv_1 = 2.3*(2*kp_1)^.5;
kp_2 = 0; % 20000,10000,  5000,  3000,2000,1000,500,200,100,0
kv_2 = 1*(2*kp_2)^.5;
kp_3 = 0; % 5500
kv_3 = 0.8*(2*kp_3)^.5;
alpha = 0.001; 
Fe = 0;

T = 0.4;
dt = 0.001;
t = (0:dt:T)';
demo_x = t;
demo_y = -(demo_x - sin(2*pi*demo_x.^3).*cos(2*pi*demo_x.^3).*exp(demo_x.^4));
x_end = demo_x(end) + 0.0;
y_end = demo_y(end) + 0.0;

xdd = zeros(length(t),1);
xd = zeros(length(t),1);
x = zeros(length(t),1);
xdd(1,1) = 0;
xd(1,1) = 0;
x(1,1) = 0.0;

ydd = zeros(length(t),1);
yd = zeros(length(t),1);
y = zeros(length(t),1);
ydd(1,1) = 0;
yd(1,1) = 0;
y(1,1) = 0.5;
R_search_start = 1;
R_search_end = length(t);%length(x)/10;
index_data = zeros(2,1);
index_data(1,1) = 1; 
index_data(2,1) = +inf; 
R_search = 100;
steps = 1500;
a = steps; %Õ÷‘≤≥§÷·
b = 1; %Õ÷‘≤∂Ã÷·

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
    
    ydd(i+1) = kp_1*y_a - kv_1*yd(i) + k_e*(kp_2*y_r - kv_2*yd(i) + kp_3*y_t + kv_3*yd(i)) + Fe; % e+r+t
    yd(i+1) = yd(i) + ydd(i)*dt;
    y(i+1) = y(i) + yd(i)*dt;
    
    xdd(i+1) = kp_1*x_a - kv_1*xd(i) + k_e*(kp_2*x_r - kv_2*xd(i) + kp_3*x_t + kv_3*xd(i)) + Fe; % e+r+t
    xd(i+1) = xd(i) + xdd(i)*dt;
    x(i+1) = x(i) + xd(i)*dt;
    
    index_data(2,1) = +inf;
    
    if  rem(i,1)==0
%         plot([x(i),x(i) + x_r], [y(i),y(i) + y_r],'r');
%         hold on;
%         plot([x(i),x(i) + x_t], [y(i),y(i) + y_t]);
%         hold on;
%         plot(x(i), k_e, 'r*');
%         hold on;
    end
end

% plot(demo_x, demo_y, 'b', 'LineWidth',1);
hold on;
% plot(x, y);
hold on;
plot(x(1,1), y(1,1), 'r*'); 
hold on;
% legend('demonstration','repro1,kp\_2=20000','repro2,kp\_2=10000','repro3,kp\_2=5000','repro4,kp\_2=3000','repro5,kp\_2=2000',...
%                        'repro6,kp\_2=1000','repro7,kp\_2=500','repro8,kp\_2=200','repro9,kp\_2=100', 'repro10,kp\_2=0');
xlabel('x','fontsize',18);
ylabel('y','fontsize',18);
set(gca,'FontSize',12);
grid on;
