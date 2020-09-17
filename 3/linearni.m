%% init
clear variables; close all; clc;

%% generisanje podataka za linearni klasifikator
N = 5000; % br odbiraka

P11 = 0.6;
M11 = [1; 1];
M12 = [6; 4];
S11 = [4, 1.1; 1.1, 2];
S12 = [2, 1.1; 1.1, 4];

P22 = 0.55;
M21 = [7; -4];
M22 = [0; -5];
S21 = [3, -0.8; -0.8, 1.5];
S22 = [5, 1; 1, 1];

X1 = zeros(N,2);
X2 = zeros(N,2);

for i=1:N
    t = rand(1,2);
    switch(t(1) < P11)
        case true
            X1(i,:) = mvnrnd(M11,S11);
        case false
            X1(i,:) = mvnrnd(M12,S12);
    end
    switch(t(2) < P22)
        case true
            X2(i,:) = mvnrnd(M21,S21);
        case false
            X2(i,:) = mvnrnd(M22,S22);
    end
end
xaxis = [min([X1(:,1); X2(:,1)],[],'All'), max([X1(:,1); X2(:,1)],[],'All')];
xaxis = [xaxis(1)-0.2*(xaxis(2)-xaxis(1)),xaxis(2)+0.2*(xaxis(2)-xaxis(1))];
yaxis = [min([X1(:,2); X2(:,2)],[],'All'), max([X1(:,2); X2(:,2)],[],'All')];
yaxis = [yaxis(1)-0.2*(yaxis(2)-yaxis(1)),yaxis(2)+0.2*(yaxis(2)-yaxis(1))];


figure(1)
plot(X1(:,1),X1(:,2), 'r.'); hold on;
plot(X2(:,1),X2(:,2), 'b.');
xlim(xaxis); ylim(yaxis);
%% III iterativna metoda
NO = round(0.7*N);
NT = N - NO;
O = logical([ones(1,NO) zeros(1,NT)]);
T = logical([zeros(1,NO) ones(1,NT)]);
M1 = mean(X1(O,:));
M2 = mean(X2(O,:));
S1 = cov(X1(O,:));
S2 = cov(X2(O,:));
itt_vector = zeros(1001,3);% [s greska v0]
for s = 0:0.001:1
    V = (s*S1+(1-s)*S2)^(-1)*(M2-M1)';
    y1t = V'*X1(T,:)';
    y2t = V'*X2(T,:)';
    y1o = V'*X1(O,:)';
    y2o = V'*X2(O,:)';
    ymin = -max([y1o y2o]);
    ymax = -min([y1o y2o]);
    i = 0;
    err = zeros(1,1000);
    v0_steps = linspace(ymin, ymax, 1000);
    for v0 = v0_steps
        i = i + 1;
        err(i) = sum(~(y1o<-v0))+sum(~(y2o>-v0));
    end
    [~, i] = min(err);
    v0 = v0_steps(i);
    itt_vector(round(s*1000+1),:) = [s sum(~(y1t<-v0))+sum(~(y2t>-v0)) v0];
end
figure(2);
plot(itt_vector(:,1),itt_vector(:,2));
[err_min, i] = min(itt_vector(:,2));
s = itt_vector(i,1);
v0 = itt_vector(i,3);
V = (s*S1+(1-s)*S2)^(-1)*(M2-M1)';
x = -10:0.1:20;
y = -10:0.1:10;
h = zeros(length(x),length(y));
figure(1); hold on;
plot(xaxis,[(-V(1)*xaxis(1)-v0)/V(2),(-V(1)*xaxis(2)-v0)/V(2)],'LineWidth',2);

%% metod zeljenog izlaza

Z = [ones(1,N),-ones(1,N);X1',-X2'];
G = [ones(1,N), ones(1,N)];
W = (Z*Z')^(-1)*Z*G';
v0 = W(1); V = W(2:3);
plot(xaxis,[(-V(1)*xaxis(1)-v0)/V(2),(-V(1)*xaxis(2)-v0)/V(2)],'LineWidth',2);
G = [ones(1,N), 2*ones(1,N)];
W = (Z*Z')^(-1)*Z*G';
v0 = W(1); V = W(2:3);
plot(xaxis,[(-V(1)*xaxis(1)-v0)/V(2),(-V(1)*xaxis(2)-v0)/V(2)],'LineWidth',2);
G = [2*ones(1,N), ones(1,N)];
W = (Z*Z')^(-1)*Z*G';
v0 = W(1); V = W(2:3);
plot(xaxis,[(-V(1)*xaxis(1)-v0)/V(2),(-V(1)*xaxis(2)-v0)/V(2)],'LineWidth',2);
legend('Klasa 1','Klasa 2','Lin. klas. III it. metoda',...
       'MZI G = [1; 1]','MZI G = [2; 1]','MZI G = [1; 2]');
   
