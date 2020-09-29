%% init
clear variables; close all; clc;

%% generisanje podataka
N = 500; % br odbiraka
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

P1 = 0.5; P2 = 0.5;

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

%% kvadratna dekompozicija
Nitt = 100; 
duzine = [1 2 4 10 20]; %25 50 100 125 250 500];
iter_vektn = zeros(1,length(duzine)+1);
X = [X1; X2];
X_1 = X; index_x = 1:length(X);
X_2 = X; index_y = 1:length(X);
for i=1:2*N-1
    for j=2:2*N
        if X_1(i,1)>X_1(j,1)
            temp = X_1(i,1);
            X_1(i,1) = X_1(j,1);
            X_1(j,1) = temp;
            temp = index_x(i);
            index_x(i) = index_x(j);
            index_x(j) = temp;
        end
    end
end
for i=1:2*N-1
    for j=2:2*N
        if X_2(i,2)>X_2(j,2)
            temp = X_2(i,2);
            X_2(i,2) = X_2(j,2);
            X_2(j,2) = temp;
            temp = index_y(i);
            index_y(i) = index_y(j);
            index_y(j) = temp;
        end
    end
end

for itt = 1:Nitt
    for d = 1:length(duzine)
        lst = 1:length(X);
        Nb = 2*N/duzine(d);
        indexes = randperm(Nb);
        lst = reshape(lst,Nb,length(X)/Nb);

        K1 = lst(indexes(1:Nb/2),:);
        K1 = reshape(K1',1,numel(K1));
        K2 = lst(indexes(Nb/2+1:Nb),:);
        K2 = reshape(K2',1,numel(K2));

        
        Mk1 = mean(X(K1,:),1);Mk2 = mean(X(K2,:),1);
        Sk1 = cov(X(K1,:));Sk2 = cov(X(K2,:));
        Pk1 = P1; Pk2 = P2;
        isOver = false;
        i = 0;
        while (~isOver && i<100)
            i = i + 1;
            isOver = true;
            j = 1;
            while (j <= length(K1))
                d1= ndek(X(K1(j),:),Mk1,Sk1,Pk1);
                d2= ndek(X(K1(j),:),Mk2,Sk2,Pk2);
                dmin= min([d1,d2]);
                if (dmin == d2) && (length(K1)>1)
                    K2 = [K2 K1(j)];
                    K1 = [K1(1:j-1) K1(j+1:end)]; 
                    isOver = false; j=j-1;
                end
                j=j+1;
            end
            
            j = 1;
            while (j <= length(K2))
                d1= ndek(X(K2(j),:),Mk1,Sk1,Pk1);
                d2= ndek(X(K2(j),:),Mk2,Sk2,Pk2);
                dmin= min([d1,d2]);
                if (dmin == d1) && (length(K2)>1)
                    K1 = [K1 K2(j)];
                    K2 = [K2(1:j-1) K2(j+1:end)];
                    isOver = false; j=j-1;
                end
                j=j+1;
            end
            Mk1 = mean(X(K1,:),1);Mk2 = mean(X(K2,:),1);
            Sk1 = cov(X(K1,:));Sk2 = cov(X(K2,:));
            Pk1 = length(K1)/(N*2); Pk2 = length(K2)/(N*2);
        end
           iter_vektn(d) = iter_vektn(d)+i/Nitt;
    end
end    

function v = ndek(X,M,S,P)
    v = 0.5*(X-M)*(S^-1)*(X-M)'+0.5*log(det(S))-0.5*log(P);
end