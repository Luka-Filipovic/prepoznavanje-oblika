%% init
clear variables; close all; clc;

%% generisanje podataka za kvadratni klasifikator
N = 5000; % br odbiraka
P11 = 0.6;
M11 = [1; 1];
M12 = [9; 2];
S11 = [8, 1.1; 1.1, 1];
S12 = [3, -2.1; -2.1, 4];

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

%% kvadratni klasifikator
G = [ones(1,N), ones(1,N)];
Z = [ones(1,N) -ones(1,N); 
     X1',-X2'; 
     X1(:,1)'.^2 -X2(:,1)'.^2;
     2.*X1(:,1)'.*X1(:,2)' -2.*X2(:,1)'.*X2(:,2)'; 
     X1(:,2)'.^2 -X2(:,2)'.^2;];
 
W = (Z*Z')^(-1)*Z*G';
v0 = W(1);
V = W(2:3); Q = [W(4) W(5);W(5) W(6)]; 
syms xp yp;
xp = solve(v0+V(1)*xp+V(2).*yp+xp^2*Q(1,1)+xp.*yp.*2*Q(1,2)+yp.^2*Q(2,2),xp);

yp = yaxis(1):0.01:yaxis(2);
xp = eval(xp);
xplot = [xp(1,:) fliplr((xp(2,:)))];
yplot = [yp fliplr(yp)];
xplot1 = xplot(imag(xplot)==0);
yplot1 = yplot(imag(xplot(1,:))==0);
figure(1)
hold on
plot(xplot1,yplot1,'g','LineWidth',2);
legend('Klasa 1','Klasa 2','Granica kvadratnog klasifikatora');

xlim(xaxis); ylim(yaxis);