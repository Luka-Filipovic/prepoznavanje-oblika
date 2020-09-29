%% Init
clear variables; close all; clc;

%% Generisanje podataka
N = 500;

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

%% Wald
eps1 = 10^-4;
eps2 = 10^-4;
A = (1-eps1)/eps2; B = eps1/(1-eps2);
a = -log(A);
b = -log(B);
N_itt = 100;
itt = zeros(2,N_itt);
for i = 1:N_itt
    s = [0];
    isOver = false;
    while(~isOver)
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
            isOver = true;
            figure(2)
            hold on;
            plot(0:itt(1,i), s, 'r');
        end
    end
    s = [0];
    isOver = false;
    while(~isOver)
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
            isOver = true;
            figure(2)
            hold on;
            plot(0:itt(2,i), s, 'b');
        end
    end
end

plot([0 max(itt,[],'All')],[1 1].*a,'k--')
plot([0 max(itt,[],'All')],[1 1].*b,'k--')
xlim([0, max(itt,[],'All')]);