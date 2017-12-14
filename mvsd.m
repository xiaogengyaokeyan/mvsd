%% MVSD Arbitrary-Start Profile 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kp_1 = 500; 
% kv_1 = 2.3*(2*kp_1)^.5; %0.08*(2*kp_1)^.5, 0.8*(2*kp_1)^.5;
% kp_2 = 20000,10000,  5000,  3000,2000,1000,500,200,100,0
% kv_2 = 1*(2*kp_2)^.5;
% kp_3 = 0; 
% kv_3 = 0.8*(2*kp_3)^.5;
% alpha = 0.001;
% Fe = 0;
% T = 0.4;
% dt = 0.001;
% x(1,1) = 0,0.1,0.2; profile 0.0
% y(1,1) = 0.4,0.3,0.2,0.1,0; profile 0.5
% demo_y(end)-y(i)-0.0
% demo_x(end)-x(i)-0.0
% R_search = 100;
% steps = 1500;

%% MVSD Arbitrary-Attractor Profile 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kp_1 = 500; 
% kv_1 = 2.3*(2*kp_1)^.5; %0.08*(2*kp_1)^.5, 0.8*(2*kp_1)^.5;
% kp_2 = % 200000,100000,50000,  20000,  10000,5000,3000,1000,0
% kv_2 = 1*(2*kp_2)^.5;
% kp_3 = 0; 
% kv_3 = 0.8*(2*kp_3)^.5;
% alpha = 0.01; % 0.015;
% Fe = 0;
% T = 0.4;
% dt = 0.001;
% x(1,1) = 0;
% y(1,1) = 0;
% x_end = demo_x(end) - 0.1; % 0.1, 0.0, -0.1
% y_end = demo_y(end) - 0.0; % -0.1, -0.05, 0.0, 0.05, 0.1
% R_search = 100;
% steps = 2000;

%% MVSD Multi-Demostrations Profile 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kp_1 = 1000; 
% kv_1 = 2.3*(2*kp_1)^.5; %0.08*(2*kp_1)^.5, 0.8*(2*kp_1)^.5;
% kp_2 = 20000,10000,5000,3000,1000
% kv_2 = 1*(2*kp_2)^.5;
% kp_3 = 0; 
% kv_3 = 0.8*(2*kp_3)^.5;
% alpha = 0.005,0.01;
% Fe = 0;
% T = 0.4;
% dt = 0.001;
% x(1,1) = 0;
% y(1,1) = 0;
% x_end = demo_x(end) - 0.1;
% y_end = demo_y(end) - 0.01,0.02,0.03,0.04,0.05;
% R_search = 100;
% steps = 1000;

%% MVSD 0.8s trajectory track 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kp_1 = 200;
% kv_1 = 2.3*(2*kp_1)^.5; 
% kp_2 = 10000;
% kv_2 = 1*(2*kp_2)^.5;
% kp_3 = 550;
% kv_3 = 0.8*(2*kp_3)^.5;
% alpha = 0.0005; 
% Fe = 0;
% T = 0.8;
% dt = 0.001;
% x(1,1) = 0.0;
% y(1,1) = 0.0;
% x_end = demo_x(end) + 0.0;
% y_end = demo_y(end) - 0.0;
% R_search = 100;
% steps = 2000;

%% varying stiffness coefficient functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k_e = exp(-alpha*i);
% k_e = (i/steps-1)^2;
% k_e = (-1/steps)*i + 1;
% k_e = -(1/(steps)*i)^2 + 1;
% k_e = b/a*(a^2-i^2)^.5;
% k_ey = exp(-alpha*1/(y_a+0.000001));
% k_ex = exp(-alpha*1/(x_a+0.000001));

%%
clear; close all; clc;
kp_1 = 200;
kv_1 = 2.3*(2*kp_1)^.5;
kp_2 = 10000;
kv_2 = 1*(2*kp_2)^.5;
kp_3 = 550; 
kv_3 = 0.8*(2*kp_3)^.5;
alpha = 0.0005; 
Fe = 0;

%%
T = 0.8;
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
y(1,1) = 0.0;

R_search_start = 1;
R_search_end = length(t);
index_data = zeros(2,1);
index_data(1,1) = 1; 
index_data(2,1) = +inf; 
R_search = 100;
steps = 2000;

%%
for i = 1:steps %length(t)-1
    for j = R_search_start:R_search_end %length(x)
        s = (demo_x(j) - x(i))^2 + (demo_y(j) - y(i))^2;
        if s <= index_data(2,1)
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
%         plot([x(i),x(i) + x_r], [y(i),y(i) + y_r],'r');
%         hold on;
%         plot([x(i),x(i) + x_t/100], [y(i),y(i) + y_t/100]);
%         hold on;
%         plot(x(i), k_e, 'r*');
%         hold on;
    end
end


%%
% subplot(1,2,1);
% set(gca, 'Units', 'normalized', 'Position', [0.1 0.3 0.4 0.5])
% set(gca, 'Units', 'normalized', 'Position', [0.6 0.3 0.4 0.5])
plot(demo_x, demo_y, 'b', 'LineWidth',1);
hold on;
plot(x, y, 'r', 'LineWidth',1);
hold on;
% x1 = [0.0 0.0 0.0 0.0 0.0 0.1 0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.2];
% y1 = [0.0 0.1 0.2 0.3 0.4 0.0 0.1 0.2 0.3 0.4 0.0 0.1 0.2 0.3 0.4];
% x1 = [demo_x(end)+0.1 demo_x(end)+0.1 demo_x(end)+0.1 demo_x(end)+0.1 demo_x(end)+0.1...
%       demo_x(end)+0.0 demo_x(end)+0.0 demo_x(end)+0.0 demo_x(end)+0.0 demo_x(end)+0.0...
%       demo_x(end)-0.1 demo_x(end)-0.1 demo_x(end)-0.1 demo_x(end)-0.1 demo_x(end)-0.1];
% y1 = [demo_y(end)-0.1 demo_y(end)-0.05 demo_y(end)+0.0 demo_y(end)+0.05 demo_y(end)+0.1...
%       demo_y(end)-0.1 demo_y(end)-0.05 demo_y(end)+0.0 demo_y(end)+0.05 demo_y(end)+0.1...
%       demo_y(end)-0.1 demo_y(end)-0.05 demo_y(end)+0.0 demo_y(end)+0.05 demo_y(end)+0.1];
% plot(x1, y1, 'r*');
% plot(x(1,1), y(1,1), 'r*'); 
% plot(x_end, y_end, 'r*');
hold on;
% legend('demonstration','repro1,kp\_2=200000','repro2,kp\_2=100000','repro3,kp\_2=50000','repro4,kp\_2=20000','repro5,kp\_2=10000',...
%                        'repro6,kp\_2=5000','repro7,kp\_2=3000','repro8,kp\_2=1000','repro9,kp\_2=0');
% legend('demonstration','repro1,(0.1,-0.1)','repro2,(0.1,-0.05)','repro3,(0.1,0.0)','repro4,(0.1,0.05)','repro5,(0.1,0.1)',...
%                        'repro6,(0.0,0.1)','repro7,(0.0,0.05)','repro8,(0.0,0.0)','repro9,(0.0,-0.05)','repro10,(0.0,-0.1)',...
%                        'repro11,(-0.1,-0.1)','repro12,(-0.1,-0.05)','repro13,(-0.1,0.0)','repro14,(-0.1,0.05)','repro15,(-0.1,0.1)');
% title('Reproduction of Demonstration with MVSD');
title('Reproduction of 2-D Demonstration with MVSD');
xlabel('x','fontsize',18);
ylabel('y','fontsize',18);
% xlim([0,0.4]);
% set(gca,'XTick',0:0.1:0.5);
% set(gca,'YTick',0:200:1000);
set(gca,'FontSize',12);
% axis equal;
grid on;
set(gca, 'GridLineStyle' ,'--');
