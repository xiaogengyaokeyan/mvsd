clear; close all; clc;
%% demox, demoy, demoz parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kp_1 = 500; % 200
% kv_1 = 2.3*(2*kp_1)^.5;
% kp_2 = 40000; % 40000
% kv_2 = 1*(2*kp_2)^.5;
% kp_3 = 0; % 50
% kv_3 = 0.8*(2*kp_3)^.5;
% alpha = 0; % 0.0001
% Fe = 0;
% R_search = 1000;
% steps = 2000;

%% demox1, demoy1, demoz1 parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kp_1 = 200; % 200
% kv_1 = 2.3*(2*kp_1)^.5;
% kp_2 = 4900; % 4000
% kv_2 = 1*(2*kp_2)^.5;
% kp_3 = 1; % 1
% kv_3 = 0.8*(2*kp_3)^.5;
% alpha = 0; % 0.0001
% Fe = 0;
% R_search = 1600;
% steps = 2000;

%% demox2, demoy2, demoz2 parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kp_1 = 400; % 200
% kv_1 = 2.3*(2*kp_1)^.5;
% kp_2 = 20000; % 4000
% kv_2 = 1*(2*kp_2)^.5;
% kp_3 = 40; % 1
% kv_3 = 0.8*(2*kp_3)^.5;
% alpha = 0; % 0.0001
% Fe = 0;
% R_search = 1000;
% steps = 1000;

%% demox3, demoy3, demoz3 parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kp_1 = 500; % 200
% kv_1 = 2.3*(2*kp_1)^.5;
% kp_2 = 10000; % 4000
% kv_2 = 1*(2*kp_2)^.5;
% kp_3 = 0; % 1
% kv_3 = 0.8*(2*kp_3)^.5;
% alpha = 0; % 0.0001
% Fe = 0;
% R_search = 500;
% steps = 1500;

%% demox4, demoy4, demoz4 parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kp_1 = 50; 
% kv_1 = 2.3*(2*kp_1)^.5;
% kp_2 = 10000; 
% kv_2 = 1*(2*kp_2)^.5;
% kp_3 = 200; 
% kv_3 = 0.8*(2*kp_3)^.5;
% alpha = 0;
% Fe = 0;
% R_search = 200;
% steps = 1190;

%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp_1 = 50; % 100,50
kv_1 = 2.3*(2*kp_1)^.5;
kp_2 = 10000; % 10000
kv_2 = 1*(2*kp_2)^.5;
kp_3 = 120; %200,120
kv_3 = 1*(2*kp_3)^.5;
alpha = 0;
Fe = 0;
k = 0.000248; % 0.001,0.000235

%% initiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 0.4;
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

R_search_start = 1;
R_search_end = length(t);
index_data = zeros(2,1);
index_data(1,1) = 1; 
index_data(2,1) = +inf; 
R_search = 800;% 200, 800
steps = 11866;% 1190 4700

x_end = demo_x(end) + 0.0;
y_end = demo_y(end) + 0.0;
z_end = demo_z(end) + 0.0;

xdd = zeros(steps,1);
xd = zeros(steps,1);
x = zeros(steps,1);
xdd(1,1) = 0;
xd(1,1) = 0;
x(1,1) = demo_x(1,1);

ydd = zeros(steps,1);
yd = zeros(steps,1);
y = zeros(steps,1);
ydd(1,1) = 0;
yd(1,1) = 0;
y(1,1) = demo_y(1,1);

zdd = zeros(steps,1);
zd = zeros(steps,1);
z = zeros(steps,1);
zdd(1,1) = 0;
zd(1,1) = 0;
z(1,1) = demo_z(1,1);

x_r = zeros(steps,1);
y_r = zeros(steps,1);
z_r = zeros(steps,1);

x_t = zeros(steps,1);
y_t = zeros(steps,1);
z_t = zeros(steps,1);

x_a = zeros(steps,1);
y_a = zeros(steps,1);
z_a = zeros(steps,1);

k_e = zeros(steps,1);
%% search and integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:steps %length(t)-1
    for j = R_search_start:R_search_end %length(x)
        s = (demo_x(j) - x(i))^2 + (demo_y(j) - y(i))^2 + (demo_z(j) - z(i))^2;
        if s <= index_data(2,1) % '='forcing data moving ahead
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
    
    index_data(1,1)
    x_r(i) = demo_x(index_data(1,1))-x(i);
    y_r(i) = demo_y(index_data(1,1))-y(i);
    z_r(i) = demo_z(index_data(1,1))-z(i);
    
    if index_data(1,1) < length(t)
        x_t(i) = demo_x(index_data(1,1)+1) - demo_x(index_data(1,1));
        y_t(i) = demo_y(index_data(1,1)+1) - demo_y(index_data(1,1));
        z_t(i) = demo_z(index_data(1,1)+1) - demo_z(index_data(1,1));
        t_mod = (y_t(i)^2+x_t(i)^2+z_t(i)^2)^.5;
        x_t(i) = x_t(i)/t_mod;
        y_t(i) = y_t(i)/t_mod;
        z_t(i) = z_t(i)/t_mod;
    end
    
    x_a(i) = x_end - x(i);
    y_a(i) = y_end - y(i);
    z_a(i) = z_end - z(i);
    k_e(i) = exp(-alpha*i);
    
    %xdd(i+1) = kp_1*x_a - kv_1*xd(i) + k_ex*(kp_2*(demo_x(i)-x(i)) - kv_2*xd(i)) + Fe; % e+i
    %xdd(i+1) = kp_1*x_t - kv_1*xd(i) + k_ex*(kp_2*x_r - kv_2*xd(i)) + Fe; % t+r, ok
    xdd(i+1) = k*(kp_1*x_a(i) - kv_1*xd(i) + k_e(i)*(kp_2*x_r(i) - kv_2*xd(i) + kp_3*x_t(i) - kv_3*xd(i)) + Fe); % e+r+t
    %xdd(i+1) = kp_1*(demo_x(i)-x(i)) - kv_1*xd(i) + Fe; % vsd
    xd(i+1) = xd(i) + xdd(i)*dt;
    x(i+1) = x(i) + xd(i)*dt;
    
    %ydd(i+1) = kp_1*y_a - kv_1*yd(i) + k_ey*(kp_2*(demo_y(i)-y(i)) - kv_2*yd(i)) + Fe; % e+i
    %ydd(i+1) = kp_1*y_t - kv_1*yd(i) + k_ey*(kp_2*y_r - kv_2*yd(i)) + Fe; % t+r, ok
    ydd(i+1) = k*(kp_1*y_a(i) - kv_1*yd(i) + k_e(i)*(kp_2*y_r(i) - kv_2*yd(i) + kp_3*y_t(i) - kv_3*yd(i)) + Fe); % e+r+t
    %ydd(i+1) = kp_1*(demo_y(i)-y(i)) - kv_1*yd(i) + Fe; % vsd
    yd(i+1) = yd(i) + ydd(i)*dt;
    y(i+1) = y(i) + yd(i)*dt;
    
    zdd(i+1) = k*(kp_1*z_a(i) - kv_1*zd(i) + k_e(i)*(kp_2*z_r(i) - kv_2*zd(i) + kp_3*z_t(i) - kv_3*zd(i)) + Fe); % e+r+t
    zd(i+1) = zd(i) + zdd(i)*dt;
    z(i+1) = z(i) + zd(i)*dt;
    
    index_data(2,1) = +inf;
    
    if  rem(i,10)==0
%         plot3([x(i),x(i) + x_r(i)], [y(i),y(i) + y_r(i)], [z(i),z(i) + z_r(i)],'r');
%         hold on;
%         plot3([x(i),x(i) + x_t(i)/50], [y(i),y(i) + y_t(i)/50], [z(i),z(i) + z_t(i)/50]);
%         hold on;
%         plot(x(i), k_e(i), 'r*');
%         hold on;
    end
end

%% plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
plot3(demo_x, demo_y,demo_z, 'b', 'LineWidth',1);
hold on;
plot3(x, y, z, 'r', 'LineWidth',1);
hold on;
plot3(demo_x(1,1), demo_y(1,1), demo_z(1,1), 'b*'); 
plot3(demo_x(end), demo_y(end), demo_z(end), 'b*');
plot3(x(1,1), y(1,1), z(1,1), 'r*'); 
plot3(x_end, y_end, z_end, 'r*');
hold on;
title('Reproduction of Demonstration with MVSD\_3d');
xlabel('x','fontsize',18);
ylabel('y','fontsize',18);
zlabel('z','fontsize',18);
set(gca,'FontSize',12);
grid on;
% axis equal;

figure;
plot(t(end)/(steps+1):t(end)/(steps+1):t(end),x, 'r');
hold on;
plot(t(end)/(steps+1):t(end)/(steps+1):t(end),y, 'g');
hold on;
plot(t(end)/(steps+1):t(end)/(steps+1):t(end),z, 'b');
hold on;
plot(t,demo_x);
hold on;
plot(t,demo_y);
hold on;
plot(t,demo_z);
hold on;
title('X, Y, Z Dimension Reproduction of Demonstration with MVSD\_3d');
xlabel('t','fontsize',18);
ylabel('x, y, z','fontsize',18);
set(gca,'FontSize',12);
grid on;
% axis equal;

figure;
plot(t(end)/(steps):t(end)/(steps):t(end),x_r, 'r');
hold on;
plot(t(end)/(steps):t(end)/(steps):t(end),y_r, 'g');
hold on;
plot(t(end)/(steps):t(end)/(steps):t(end),z_r, 'b');
hold on;
title('X, Y, Z Dimension Reproduction Error');
xlabel('t','fontsize',18);
ylabel('error','fontsize',18);
set(gca,'FontSize',12);
grid on;
% axis equal;

figure;
r = (x_r.^2 + y_r.^2 + z_r.^2).^.5;
plot(t(end)/(steps):t(end)/(steps):t(end),r,'r');
title('Reproduction Error');
xlabel('t','fontsize',18);
ylabel('error','fontsize',18);
set(gca,'FontSize',12);
grid on;
% axis equal;

figure;
plot(t(end)/(steps+1):t(end)/(steps+1):t(end),xdd, 'r');

sum(r)/steps
