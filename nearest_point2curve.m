
clear all; clc;

Xe = 1.1;%0.1;
Ye = 1.105;%0.1;
Ze = 1.102;%0.6;
index_data = zeros(2,1);
index_data(1,1) = 1; 
index_data(2,1) = +inf; 
R_search = 1000;

%T = 7;
%Xt = (0:.001:T)';
%Yt = Xt - sin(2*pi*Xt.^3).*cos(2*pi*Xt.^3).*exp(Xt.^4);

load('Data/x_1.mat');
load('Data/x_2.mat');
load('Data/x_3.mat');

Xt = x_1;
Yt = x_2;
Zt = x_3;

R_search_start = 1;
R_search_end = 2000;%length(Xt)/10;

t1=clock; 

for j = 1:1
for i = R_search_start:R_search_end %length(Xt)
    s = (Xt(i) - Xe)^2 + (Yt(i) - Ye)^2 + (Zt(i) - Ze)^2;
    if s < index_data(2,1)
        index_data(2,1) = s;
        index_data(1,1) = i;
    end
    R_search_start = index_data(1,1) - R_search;
    if R_search_start < 0
        R_search_start = 0;
    end
    R_search_end = index_data(1,1) + R_search;
    if R_search_end > length(Xt)
        R_search_end = length(Xt);
    end
end
Xe = Xe + 0.001;
Ye = Xe + 0.005;
Ze = Xe + 0.002;
end

t2=clock;
etime(t2,t1)

figure;
Xs = Xt(index_data(1,1));
Ys = Yt(index_data(1,1));
Zs = Zt(index_data(1,1));
disp(['shortest distance is: ',num2str(index_data(2,1))]);
disp(['query point is: ',num2str([Xe, Ye, Ze])]);
disp(['shortest point in trajectory is: ',num2str([Xs, Ys, Zs])]);
plot3(Xt, Yt, Zt,'b');
hold on;
plot3(Xe,Ye,Ze,'r*',Xs,Ys,Zs,'k*');
x = Xs:(Xe-Xs)/100:Xe;
y = Ys:(Ye-Ys)/100:Ye;
z = Zs:(Ze-Zs)/100:Ze;
plot3(x,y,z,'r--');
%axis equal;
