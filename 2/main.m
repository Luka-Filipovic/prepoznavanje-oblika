clear variables; close all; clc;

%% generisanje podataka
N = 500; % br odbiraka

P1 = 0.6;
M11 = [1; 1];
M12 = [6; 4];
S11 = [4, 1.1; 1.1, 2];
S12 = [3, -0.8; -0.8, 1.5];

P2 = 0.55;
M21 = [7; -4];
M22 = [6; 0];
S21 = [2, 1.1; 1.1, 4];
S22 = [3, 0.8; 0.8, 0.5];

X1 = zeros(2,N);
X2 = zeros(2,N);

for i=1:N
    t = rand(1,2);
    switch(t(1) < P1)
        case true
            X1(:,i) = mvnrnd(M11,S11);
        case false
            X1(:,i) = mvnrnd(M12,S12);
    end
    switch(t(2) < P2)
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

for i = 1:length(x)
    for j = 1:length(y)
        X = [x(i);y(j)];
        f11 = 1/(2*pi*det(S11)^0.5)*exp(-0.5*(X-M11)'*S11^(-1)*(X-M11));
        f12 = 1/(2*pi*det(S12)^0.5)*exp(-0.5*(X-M12)'*S12^(-1)*(X-M12));
        f1(i,j) = P1*f11 + (1-P1)*f12;
        f21 = 1/(2*pi*det(S21)^0.5)*exp(-0.5*(X-M21)'*S21^(-1)*(X-M21));
        f22 = 1/(2*pi*det(S22)^0.5)*exp(-0.5*(X-M22)'*S22^(-1)*(X-M22));
        f2(i,j) = P2*f21 + (1-P2)*f22;
        h(i,j) = log(f2(i,j))-log(f1(i,j));
    end
end

f1max = max(f1,[],'All');
f2max = max(f2,[],'All');

dk = [1 3 6 9];
for d = dk
    prag=f1max*exp(-0.5*dk);
    contour(x,y,f1',prag,'m');
    contour(x,y,f2',prag,'m');
end
hold off;

figure(2)
plot(X1(1,:),X1(2,:), 'r.'); hold on;
plot(X2(1,:),X2(2,:), 'b.');
contour(x,y,h',[0 0]);
hold off;
x_step = x(2)-x(1);
y_step = y(2)-y(1);
Eps1 =sum(sum(f1(h>=0)))*x_step*y_step;
Eps2 =sum(sum(f2(h<0)))*x_step*y_step;