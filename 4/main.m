%% init
clear variables; close all; clc;

%% generisanje podataka
N = 500; % br odbiraka

M1 = [-5; -4];
M2 = [5; 4];
M3 = [3; -2];
M4 = [-4; 3];
S1 = [4, 1.1; 1.1, 2];
S2 = [3, 1.1; 1.1, 3];
S3 = [4, -1.8; -1.8, 2];
S4 = [2, 1.1; 1.1, 4];

X1 = zeros(N,2);
X2 = zeros(N,2);
X3 = zeros(N,2);
X4 = zeros(N,2);

for i=1:N
    X1(i,:) = mvnrnd(M1,S1);
    X2(i,:) = mvnrnd(M2,S2);
    X3(i,:) = mvnrnd(M3,S3);
    X4(i,:) = mvnrnd(M4,S4);
end
xaxis = [min([X1(:,1); X2(:,1); X3(:,1); X4(:,1)],[],'All'),...
         max([X1(:,1); X2(:,1); X3(:,1); X4(:,1)],[],'All')];
xaxis = [xaxis(1)-0.2*(xaxis(2)-xaxis(1)),xaxis(2)+0.2*(xaxis(2)-xaxis(1))];
yaxis = [min([X1(:,2); X2(:,2); X3(:,2); X4(:,2)],[],'All'),...
         max([X1(:,2); X2(:,2); X3(:,2); X4(:,2)],[],'All')];
yaxis = [yaxis(1)-0.2*(yaxis(2)-yaxis(1)),yaxis(2)+0.2*(yaxis(2)-yaxis(1))];

figure(1)
plot(X1(:,1),X1(:,2), 'r.'); hold on;
plot(X2(:,1),X2(:,2), 'b.');
plot(X3(:,1),X3(:,2), 'g.');
plot(X4(:,1),X4(:,2), 'c.');
xlim(xaxis); ylim(yaxis);