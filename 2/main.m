%% Init
clear variables; close all; clc;

%% Generisanje podataka
N = 500; % br odbiraka

P11 = 0.6;
M11 = [1; 1];
M12 = [6; 4];
S11 = [4, 1.1; 1.1, 2];
S12 = [3, -0.8; -0.8, 1.5];

P22 = 0.55;
M21 = [7; -4];
M22 = [6; 0];
S21 = [2, 1.1; 1.1, 4];
S22 = [3, 0.8; 0.8, 0.5];



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

for i = 1:length(x)
    for j = 1:length(y)
        X = [x(i);y(j)];
        f11 = 1/(2*pi*det(S11)^0.5)*exp(-0.5*(X-M11)'*S11^(-1)*(X-M11));
        f12 = 1/(2*pi*det(S12)^0.5)*exp(-0.5*(X-M12)'*S12^(-1)*(X-M12));
        f1(i,j) = P11*f11 + (1-P11)*f12;
        f21 = 1/(2*pi*det(S21)^0.5)*exp(-0.5*(X-M21)'*S21^(-1)*(X-M21));
        f22 = 1/(2*pi*det(S22)^0.5)*exp(-0.5*(X-M22)'*S22^(-1)*(X-M22));
        f2(i,j) = P22*f21 + (1-P22)*f22;
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
%% Klasifikator minimalne greske
figure(2)
plot(X1(1,:),X1(2,:), 'r.'); hold on;
plot(X2(1,:),X2(2,:), 'b.');
contour(x,y,h',[0 0],'DisplayName','Klasifikator minimalne greske');
x_step = x(2)-x(1);
y_step = y(2)-y(1);
Eps1 =sum(f1(h>=0),'All')*x_step*y_step;
Eps2 =sum(f2(h<0),'All')*x_step*y_step;

Mx = zeros(2,2);
for i=1:N
    [~,ix] = min(abs(x-X1(1,i)));
    [~,iy] = min(abs(y-X1(2,i)));
    if (h(ix,iy) <= 0)
        Mx(1,1) = Mx(1,1) + 1;
    else
        Mx(1,2) = Mx(1,2) + 1;
    end
    [~,ix] = min(abs(x-X2(1,i)));
    [~,iy] = min(abs(y-X2(2,i)));
    if (h(ix,iy) >= 0)
        Mx(2,2) = Mx(2,2) + 1;
    else
        Mx(2,1) = Mx(2,1) + 1;
    end
end
T = array2table(Mx,'VariableNames',{'ω1_p','ω2_p'},'RowName',...
    {'ω1_t  |','ω2_t  |'}); 
disp('Matrica konfuzije:');
disp(T);

disp('Dobijene greske:');
disp([Mx(1,2)/N, Mx(2,1)/N]);
disp('Teorijske greske:');
disp([Eps1,Eps2]);
%% Klasifikator minimalne cene
c11 = 0; c22 = 0;
c12 = 1; c21 = 10*c12;

th = log((c21-c11)/(c12-c22));

contour(x,y,h',[th th],'g','DisplayName','Klasifikator minimalne cene');
hold off;
x_step = x(2)-x(1);
y_step = y(2)-y(1);
Eps1_kmc =sum(f1(h>=th),'All')*x_step*y_step;
Eps2_kmc =sum(f2(h<th),'All')*x_step*y_step;
legend;

Mx = zeros(2,2);
for i=1:N
    [~,ix] = min(abs(x-X1(1,i)));
    [~,iy] = min(abs(y-X1(2,i)));
    if (h(ix,iy) <= th)
        Mx(1,1) = Mx(1,1) + 1;
    else
        Mx(1,2) = Mx(1,2) + 1;
    end
    [~,ix] = min(abs(x-X2(1,i)));
    [~,iy] = min(abs(y-X2(2,i)));
    if (h(ix,iy) >= th)
        Mx(2,2) = Mx(2,2) + 1;
    else
        Mx(2,1) = Mx(2,1) + 1;
    end
end
T = array2table(Mx,'VariableNames',{'ω1_p','ω2_p'},'RowName',...
    {'ω1_t  |','ω2_t  |'}); 
disp('Matrica konfuzije KMC:');
disp(T);

disp('Dobijene greske:');
disp([Mx(1,2)/N, Mx(2,1)/N]);
disp('Teorijske greske:');
disp([Eps1_kmc,Eps2_kmc]);

%% Wald

P11 = 0.6;
M11 = [1; 1];
M12 = [6; 4];
S11 = [4, 1.1; 1.1, 2];
S12 = [3, -0.8; -0.8, 1.5];

P22 = 0.55;
M21 = [1.5; 0.5];
M22 = [7; 3];
S21 = [2, 1.1; 1.1, 4];
S22 = [3, 0.8; 0.8, 0.5];


eps1 = 10^-4;
eps2 = 10^-4;
A = (1-eps1)/eps2; B = eps1/(1-eps2);
a = -log(A);
b = -log(B);
N_itt = 100;
itt = zeros(2,N_itt);
for i = 1:N_itt
    s = [0];
    wald = false;
    while(~wald)
        t = rand(1,1);
        switch(t < P11)
            case true
               X = mvnrnd(M11,S11)';
            case false
               X = mvnrnd(M12,S12)';
        end
        f11 = 1/(2*pi*det(S11)^0.5)*exp(-0.5*(X-M11)'*S11^(-1)*(X-M11));
        f12 = 1/(2*pi*det(S12)^0.5)*exp(-0.5*(X-M12)'*S12^(-1)*(X-M12));
        f1 = P11*f11 + (1-P11)*f12;

        f21 = 1/(2*pi*det(S21)^0.5)*exp(-0.5*(X-M21)'*S21^(-1)*(X-M21));
        f22 = 1/(2*pi*det(S22)^0.5)*exp(-0.5*(X-M22)'*S22^(-1)*(X-M22));
        f2 = P22*f21 + (1-P22)*f22;
        
        h = log(f2)-log(f1);

        itt(1,i) = itt(1,i) + 1;
        s = [s, s(end) + h];
        if (s(end) < a || s(end) > b)
            wald = true;
            figure(4)
            hold on;
            plot(0:itt(1,i), s, 'r');
        end
    end
    s = [0];
    wald = false;
    while(~wald)
        t = rand(1,1);
        switch(t < P22)
            case true
                X = mvnrnd(M21,S21)';
            case false
                X = mvnrnd(M22,S22)';
        end
        f11 = 1/(2*pi*det(S11)^0.5)*exp(-0.5*(X-M11)'*S11^(-1)*(X-M11));
        f12 = 1/(2*pi*det(S12)^0.5)*exp(-0.5*(X-M12)'*S12^(-1)*(X-M12));
        f1 = P11*f11 + (1-P11)*f12;

        f21 = 1/(2*pi*det(S21)^0.5)*exp(-0.5*(X-M21)'*S21^(-1)*(X-M21));
        f22 = 1/(2*pi*det(S22)^0.5)*exp(-0.5*(X-M22)'*S22^(-1)*(X-M22));
        f2 = P22*f21 + (1-P22)*f22;
        
        h = log(f2)-log(f1);

        itt(2,i) = itt(2,i) + 1;
        s = [s, s(end) + h];
        if (s(end) < a || s(end) > b)
            wald = true;
            figure(4)
            hold on;
            plot(0:itt(2,i), s, 'b');
        end
    end
end

plot([0 max(itt,[],'All')],[1 1].*a,'k--')
plot([0 max(itt,[],'All')],[1 1].*b,'k--')
xlim([0, max(itt,[],'All')]);