function mvsd_dmp
% Dynamical movement primitive (DMP) with locally weighted regression (LWR).
% @ Author:gengpeng
% @ Date:2017/9/30
% @ For:Multi-virtual-spring-damper

%% DMP Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dmp.nbWeights = 8;% Number of activation functions (gaussian)
dmp.kp = 50;% Stiffness gain
dmp.kv = (2*dmp.kp)^.5;% Damping gain (with ideal underdamped damping ratio)
dmp.alpha = 1.0;% Decay factor

%% Loading Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Data/Data.mat');
% Reproduction parameters
repro.g = Data{1,1}.x(:,end);
repro.x0 = Data{1,1}.x(:,1);
repro.xd0 = Data{1,1}.xd(:,1);
repro.xdd0 = Data{1,1}.xdd(:,1);
repro.dt = Data{1,1}.dt;
repro.nbData = length(Data{1,1}.t(1,:));

nbSamples = 4;% Number of demonstrations
nbVar = length(Data{1,1}.x(:,1));% Number of dimensions, motion variables [x1,x2,x3] 
dt = Data{1,1}.dt;% Duration of time step
for i = 1:nbSamples
    dmp.g(:,i) = Data{i,1}.x(:,end);% Target of the DMP
    dmp.x0(:,i) = Data{i,1}.x(:,1);% Start of the DMP
    nbData(i) = length(Data{i,1}.t(1,:));% Length of each trajectory
end
dmp.g(:,nbSamples+1) = repro.g;% Target of the DMP
dmp.x0(:,nbSamples+1) = repro.x0;% Start of the DMP
nbData(nbSamples+1) = repro.nbData;% Length of each trajectory

for i = 1:nbSamples
    demo_x{i,1} = Data{i,1}.x;
    demo_xd{i,1} = Data{i,1}.xd;
    demo_xdd{i,1} = Data{i,1}.xdd;
end

%% Setting the basis functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute decay term: s
for i = 1:nbSamples + 1
    s{i,1}(1,1) = 1;% Initialization of decay term
    for t = 2:nbData(i)
        s{i,1}(1,t) = s{i,1}(1,t-1) - dmp.alpha * s{i,1}(1,t-1) * dt;% Update of decay term (ds/dt=-alpha s)
    end
end

% Compute Mu(c) and sigma(h) of the weight
for i = 1:nbSamples + 1
    TimingSep = linspace(min(s{i,1}),max(s{i,1}),dmp.nbWeights+1);
    for j = 1:dmp.nbWeights
        idtmp = find(s{i,1} >= TimingSep(j) & s{i,1} < TimingSep(j+1));
        c(i,j) = mean(s{i,1}(1,idtmp));
        h(i,j) = 1.0/2E-3;
    end
end

% Compute activation weight: w
for i = 1:nbSamples + 1
    for j = 1:dmp.nbWeights
        w{i,1}(j,:) = exp(-0.5 * h(i,j) * (s{i,1} - repmat(c(i,j),1,nbData(i))) .^2);
    end
end

% Normalization of the weighting term: g
for i = 1:nbSamples + 1
    g{i,1} = w{i,1} ./ repmat(sum(w{i,1}), dmp.nbWeights, 1);
end

% Compute the shape function from demonstrations: F
for i = 1:nbSamples
    F{i,1} = (demo_xdd{i,1} - (dmp.kp * (repmat(dmp.g(:,i),1,nbData(i)) - demo_x{i,1}) - dmp.kv * demo_xd{i,1})) ./ repmat(s{i,1},nbVar,1);
end

%% Learning the parameter theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = [];
Y = [];
tmpW = [];
for i = 1:nbSamples
    X = [X;ones(nbData(i),1)];
    Y = [Y;F{i,1}'];
    tmpW = [tmpW w{i,1}];
end
for i = 1:dmp.nbWeights
    W = diag(tmpW(i,:));
    MuF(:,i) = (X'*W*X \ X'*W*Y)';
end

%% Reproduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
repro.F = (MuF * g{nbSamples+1,1});% .* repmat((dmp.g(:,nbSamples+1) - dmp.x0(:,nbSamples+1)),1,repro.nbData);
x = repro.x0;
xd = repro.xd0;
xdd = repro.xdd0;
for t = 1:repro.nbData
    repro.x(:,t) = x;
    repro.xd(:,t) = xd;
    repro.xdd(:,t) = xdd;
    
    xdd = dmp.kp * (repro.g-x) - dmp.kv*xd + repro.F(:,t)*s{nbSamples+1}(1,t);
    xd = xd + xdd*repro.dt;
    x = x + xd*repro.dt;
end

%% Save and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
save('Data/repro.mat','repro');

figure;
for i = 1:nbSamples
    plot(Data{i,1}.x(1,:),Data{i,1}.x(2,:));
    hold on;
end
hold on;
plot(repro.x(1,:),repro.x(2,:));

