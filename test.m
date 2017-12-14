clear; close all; clc;
%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp = 300; % 800
kv = 1*(2*kp)^.5; % 1*(2*kp)^.5

%% initiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.002;

load('Data/demot4');
load('Data/demox4.mat');
load('Data/demoy4.mat');
load('Data/demoz4.mat');
t = t_1;
demo_x = x_1;
demo_y = y_1;
demo_z = z_1;

R_search = 1000;% 10, 50
steps = 10000;% 40000, 10000
n_steps_forward = 80; % 10, 40
R_search_start = 1;
R_search_end = R_search; 
index_data = zeros(2,1);
index_data(1,1) = 1; 
index_data(2,1) = +inf; 

x_end = demo_x(end) + 0.0;
y_end = demo_y(end) + 0.0;
z_end = demo_z(end) + 0.0;

xdd = 0;
xd = 0;
x = demo_x(1,1);

ydd = 0;
yd = 0;
y = demo_y(1,1);

zdd = 0;
zd = 0;
z = demo_z(1,1);

x_d = demo_x(1,1);
y_d = demo_y(1,1);
z_d = demo_z(1,1);

repro_x = zeros(steps,1);
repro_y = zeros(steps,1);
repro_z = zeros(steps,1);

repro_x(1,1) = x;
repro_y(1,1) = y;
repro_z(1,1) = z;

repro_xdd = zeros(steps,1);
repro_ydd = zeros(steps,1);
repro_zdd = zeros(steps,1);

repro_xdd(1,1) = xdd;
repro_ydd(1,1) = ydd;
repro_zdd(1,1) = zdd;

%% search and integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% while 1 
for i = 1:steps
    for j = R_search_start:R_search_end 
        s = (demo_x(j) - x)^2 + (demo_y(j) - y)^2 + (demo_z(j) - z)^2;
        if s <= index_data(2,1) % '='forcing data moving ahead
            index_data(2,1) = s;
            index_data(1,1) = j;
        end
    end
    R_search_start = index_data(1,1) - R_search;
    if R_search_start < 1
        R_search_start = 1;
    end
    R_search_end = index_data(1,1) + R_search;
    if R_search_end > length(t)
        R_search_end = length(t);
    end
    
    if index_data(1,1)+n_steps_forward <= length(t)
        x_d = demo_x(index_data(1,1)+n_steps_forward);
        y_d = demo_y(index_data(1,1)+n_steps_forward);
        z_d = demo_z(index_data(1,1)+n_steps_forward);
    else
        x_d = demo_x(length(t));
        y_d = demo_y(length(t));
        z_d = demo_z(length(t));
    end
    
%     index_data(1,1)
    
   
    xdd = kp*(x_d-x) - kv*xd; % vsd_time_invariant
    xd = xd + xdd*dt;
    x = x + xd*dt;
      
    ydd = kp*(y_d-y) - kv*yd; % vsd_time_invariant
    yd = yd + ydd*dt;
    y = y + yd*dt;
    
    zdd = kp*(z_d-z) - kv*zd; % vsd_time_invariant
    zd = zd + zdd*dt;
    z = z + zd*dt;
    
    repro_x(i+1) = x;
    repro_y(i+1) = y;
    repro_z(i+1) = z;
    
    repro_xdd(i+1) = xdd;
    repro_ydd(i+1) = ydd;
    repro_zdd(i+1) = zdd;
    
    index_data(2,1) = +inf;
end

%% plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot3(demo_x, demo_y, demo_z, 'b', 'LineWidth',1);
hold on;
plot3(demo_x(1,1), demo_y(1,1), demo_z(1,1), 'b*'); 
hold on;
plot3(demo_x(end), demo_y(end), demo_z(end), 'b*');
hold on;
plot3(repro_x, repro_y, repro_z, 'r', 'LineWidth',1);
hold on;
plot3(repro_x(1,1), repro_y(1,1), repro_z(1,1), 'r*');
hold on;
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
plot(t(end)/(steps+1):t(end)/(steps+1):t(end),repro_xdd, 'r');
hold on;
plot(t(end)/(steps+1):t(end)/(steps+1):t(end),repro_ydd, 'g');
hold on;
plot(t(end)/(steps+1):t(end)/(steps+1):t(end),repro_zdd, 'b');
hold on;
