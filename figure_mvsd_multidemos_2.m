clear; close all; clc;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp_1 = 1000; 
kv_1 = 2.3*(2*kp_1)^.5; 
kp_2 = 1000; % 50000, 30000 10000, 5000, 3000, 1000 
kv_2 = 1*(2*kp_2)^.5;
kp_3 = 50000; % 1000, 3000, 5000, 10000, 30000, 50000
kv_3 = 1*(2*kp_3)^.5;
kp_4 = 0; 
kv_4 = 0.8*(2*kp_4)^.5;
kp_5 = 0; 
kv_5 = 0.8*(2*kp_5)^.5;
alpha = 0.001; % 0.005
Fe = 0;

T = 0.4;
dt = 0.001;
t = (0:dt:T)';
demo_x = t;
demo1_y = -(demo_x - sin(2*pi*demo_x.^3).*cos(2*pi*demo_x.^3).*exp(demo_x.^4));
demo2_y = -0.2 - (demo_x - sin(2*pi*demo_x.^4).*cos(2*pi*demo_x.^4).*exp(demo_x.^5));
demo3_y = -0.3*ones(length(demo_x),1);

x_end = demo_x(end) + 0.0;
y_end = (demo1_y(end)+demo3_y(end))/2 - 0.0;

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
y(1,1) = (demo1_y(1)+demo3_y(1))/2;

R1_search_start = 1;
R1_search_end = length(t);%length(x)/10;
R2_search_start = 1;
R2_search_end = length(t);
R3_search_start = 1;
R3_search_end = length(t);
index_data = zeros(6,1);
index_data(1,1) = 1; 
index_data(2,1) = +inf; 
index_data(3,1) = 1; 
index_data(4,1) = +inf;
index_data(5,1) = 1; 
index_data(6,1) = +inf;
R_search = 100;
steps = 10000;
a = steps; %Õ÷‘≤≥§÷·
b = 1; %Õ÷‘≤∂Ã÷·

for i = 1:steps %length(t)-1
    for j = R1_search_start:R1_search_end %length(x)
        s1 = (demo_x(j) - x(i))^2 + (demo1_y(j) - y(i))^2;
        s2 = (demo_x(j) - x(i))^2 + (demo2_y(j) - y(i))^2;
        s3 = (demo_x(j) - x(i))^2 + (demo3_y(j) - y(i))^2;
        if s1 < index_data(2,1)
            index_data(2,1) = s1;
            index_data(1,1) = j;
        end
        if s2 < index_data(4,1)
            index_data(4,1) = s2;
            index_data(3,1) = j;
        end
        if s3 < index_data(6,1)
            index_data(6,1) = s3;
            index_data(5,1) = j;
        end
        R1_search_start = index_data(1,1) - R_search;
        R2_search_start = index_data(3,1) - R_search;
        R3_search_start = index_data(5,1) - R_search;
        if R1_search_start < 1
            R1_search_start = 1;
        end
        if R2_search_start < 1
            R2_search_start = 1;
        end
        if R3_search_start < 1
            R3_search_start = 1;
        end
        R1_search_end = index_data(1,1) + R_search;
        R2_search_end = index_data(3,1) + R_search;
        R3_search_end = index_data(5,1) + R_search;
        if R1_search_end > length(t)
            R1_search_end = length(t);
        end
        if R2_search_end > length(t)
            R2_search_end = length(t);
        end
        if R3_search_end > length(t)
            R3_search_end = length(t);
        end
    end
    
    %index_data(1,1)
    y_r1 = demo1_y(index_data(1,1))-y(i);
    x_r1 = demo_x(index_data(1,1))-x(i);
    y_r2 = demo2_y(index_data(3,1))-y(i);
    x_r2 = demo_x(index_data(3,1))-x(i);
    y_r3 = demo3_y(index_data(5,1))-y(i);
    x_r3 = demo_x(index_data(5,1))-x(i);
    
    y_t1 = -x_r1/(y_r1+0.00001);
    x_t1 = 1;
    y_t2 = -x_r2/(y_r2+0.00001);
    x_t2 = 1;
    y_t3 = -x_r3/(y_r3+0.00001);
    x_t3 = 1;
    t1_mod = (y_t1^2+x_t1^2)^.5;
    t2_mod = (y_t2^2+x_t2^2)^.5;
    t3_mod = (y_t3^2+x_t3^2)^.5;
    y_t1 = y_t1/t1_mod;
    x_t1 = x_t1/t1_mod;
    y_t2 = y_t2/t2_mod;
    x_t2 = x_t2/t2_mod;
    y_t3 = y_t3/t3_mod;
    x_t3 = x_t3/t3_mod;
    
    y_a = y_end - y(i);
    x_a = x_end - x(i);
    
    k_e = exp(-alpha*i);
    
    ydd(i+1) = kp_1*y_a - kv_1*yd(i) + k_e*(kp_2*y_r1 - kv_2*yd(i) + kp_3*y_r3 - kv_3*yd(i) + kp_4*y_t1 + kv_4*yd(i) + kp_5*y_t3 + kv_5*yd(i)) + Fe; % e+r+t
    yd(i+1) = yd(i) + ydd(i)*dt;
    y(i+1) = y(i) + yd(i)*dt;
    
    xdd(i+1) = kp_1*x_a - kv_1*xd(i) + k_e*(kp_2*x_r1 - kv_2*xd(i) + kp_3*x_r3 - kv_3*xd(i) + kp_4*x_t1 + kv_4*xd(i) + kp_5*x_t3 + kv_5*xd(i)) + Fe; % e+r+t
    xd(i+1) = xd(i) + xdd(i)*dt;
    x(i+1) = x(i) + xd(i)*dt;
    
    index_data(2,1) = +inf;
    index_data(4,1) = +inf;
    index_data(6,1) = +inf;
    
    if  rem(i,1)==0
%         plot([x(i),x(i) + x_r], [y(i),y(i) + y_r],'r');
%         hold on;
%         plot([x(i),x(i) + x_t], [y(i),y(i) + y_t]);
%         hold on;
%         plot(x(i), k_e, 'r*');
%         hold on;
    end
end
% figure;
% plot(demo_x, demo1_y,'b','LineWidth',1.0);
hold on;
% plot(demo_x, demo3_y,'r','LineWidth',1.0);
hold on;
% plot(x, y);
hold on;
plot(x_end, y_end, 'r*');
plot(x(1,1), y(1,1), 'r*');
plot(0.0, -0.35);
% legend('demonstration\_1', 'demonstration\_3', 'repro1,kp\_2=50000,kp\_3=1000','repro2,kp\_2=30000,kp\_3=1000','repro3,kp\_2=10000,kp\_3=1000',...
%                                              'repro4,kp\_2=5000,kp\_3=1000','repro5,kp\_2=3000,kp\_3=1000','repro6,kp\_2=1000,kp\_3=1000',...
%                                              'repro7,kp\_2=1000,kp\_3=3000','repro8,kp\_2=1000,kp\_3=5000','repro9,kp\_2=1000,kp\_3=10000',...
%                                              'repro10,kp\_2=1000,kp\_3=30000','repro11,kp\_2=1000,kp\_3=50000');
% title('Reproduction of Multi-Demostrations with MVSD 2');
xlabel('x','fontsize',18);
ylabel('y','fontsize',18);
set(gca,'FontSize',12);
% axis equal;
grid on;
