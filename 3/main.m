%% init
clear variables; close all; clc;

%% generisanje podataka za linearni klasifikator
N = 5000; % br odbiraka

P11 = 0.6;
M11 = [1; 1];
M12 = [6; 4];
S11 = [4, 1.1; 1.1, 2];
S12 = [3, -0.8; -0.8, 1.5];

P22 = 0.55;
M21 = [7; -4];
M22 = [0; -5];
S21 = [2, 1.1; 1.1, 4];
S22 = [5, 1; 1, 1];

X1 = zeros(2,N);
X2 = zeros(2,N);

for i=1:N
    t = rand(1,2);
    switch(t(1) < P11)
        case true
            X1(:,i) = mvnrnd(M11,S11);
        case false
            X1(:,i) = mvnrnd(M12,S12);
    end
    switch(t(2) < P22)
        case true
            X2(:,i) = mvnrnd(M21,S21);
        case false
            X2(:,i) = mvnrnd(M22,S22);
    end
end

figure(1)
plot(X1(1,:),X1(2,:), 'r.'); hold on;
plot(X2(1,:),X2(2,:), 'b.');

x = linspace(min([X1(1,:),X2(1,:)]),max([X1(1,:),X2(1,:)]),200);
y = linspace(min([X1(2,:),X2(2,:)]),max([X1(2,:),X2(2,:)]),200);

f1 = zeros(length(x),length(y));
f2 = f1;
h = f1;

%% III iterativna metoda
NO = round(0.7*N);
NT = N - NO;
M1 = mean(X1(1:NO));
M2 = mean(X2(1:NO));
S1 = var(X1(1:NO));
S2 = var(X2(1:NO));
for s = 0:0.001:1
    V = (s*S1+(1-s)*S2)^(-1)/(M2-M1);
end