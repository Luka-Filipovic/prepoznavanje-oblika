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
plot(X1(1,:),X1(2,:), '*r'); hold on;
plot(X2(1,:),X2(2,:), '*b'); hold off;