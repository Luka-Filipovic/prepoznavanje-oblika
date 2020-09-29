%% init
clear variables; close all; clc;

%% generisanje podataka
N = 500; % br odbiraka

P1 = 0.25; P2 = 0.25; P3 = 0.25; P4 = 0.25;

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

%% C-mean algoritam klasterizacije

Nitt = 100; 
duzine = [1 2 5 10 20];
itt_vect = zeros(1,length(duzine)+1);
X = [X1; X2; X3; X4];
X_1 = X; index_x = 1:length(X);
X_2 = X; index_y = 1:length(X);
for i=1:4*N-1
    for j=2:4*N
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
for i=1:4*N-1
    for j=2:4*N
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
        Nb = 4*N/duzine(d);
        indexes = randperm(Nb);
        lst = reshape(lst,Nb,length(lst)/Nb);

        K1 = lst(indexes(1:Nb/4),:);
        K1 = reshape(K1',1,numel(K1));
        K2 = lst(indexes(Nb/4+1:2*Nb/4),:);
        K2 = reshape(K2',1,numel(K2));
        K3 = lst(indexes(2*Nb/4+1:3*Nb/4),:);
        K3 = reshape(K3',1,numel(K3));
        K4 = lst(indexes(3*Nb/4+1:end),:);
        K4 = reshape(K4',1,numel(K4));     
        
        Mk1 = mean(X(K1,:),1);Mk2 = mean(X(K2,:),1);
        Mk3 = mean(X(K3,:),1);Mk4 = mean(X(K4,:),1);
        radi = true;
        i = 0;
        while (radi && i<10000)
            i = i + 1;
            radi = false;
            j = 1;
            while (j <= length(K1))
                d1= edst(X(K1(j),:),Mk1);
                d2= edst(X(K1(j),:),Mk2);
                d3= edst(X(K1(j),:),Mk3);
                d4= edst(X(K1(j),:),Mk4);
                dmin= min([d1,d2,d3,d4]);
                if (dmin == d2) && (length(K1)>1)
                    K2 = [K2 K1(j)];
                    K1 = [K1(1:j-1) K1(j+1:end)]; 
                    radi = true; j=j-1;
                    
                elseif (dmin == d3) && (length(K1)>1)
                    K3 = [K3 K1(j)];
                    K1 = [K1(1:j-1) K1(j+1:end)];
                    radi = true; j=j-1;
                    
                elseif (dmin == d4) && (length(K1)>1)
                    K4 = [K4 K1(j)];
                    K1 = [K1(1:j-1) K1(j+1:end)];
                    radi = true; j=j-1;
                end
                j=j+1;
            end
            
            j = 1;
            while (j <= length(K2))
                d1= edst(X(K2(j),:),Mk1);
                d2= edst(X(K2(j),:),Mk2);
                d3= edst(X(K2(j),:),Mk3);
                d4= edst(X(K2(j),:),Mk4);
                dmin= min([d1,d2,d3,d4]);
                if (dmin == d1) && (length(K2)>1)
                    K1 = [K1 K2(j)];
                    K2 = [K2(1:j-1) K2(j+1:end)];
                    radi = true; j=j-1;
                elseif (dmin == d3) && (length(K2)>1)
                    K3 = [K3 K2(j)];
                    K2 = [K2(1:j-1) K2(j+1:end)];
                    radi = true; j=j-1;
                elseif (dmin == d4) && (length(K2)>1)
                    K4 = [K4 K2(j)];
                    K2 = [K2(1:j-1) K2(j+1:end)];
                    radi = true; j=j-1;
                end
                j=j+1;
            end
            
            j = 1;
            while (j <= length(K3))
                d1= edst(X(K3(j),:),Mk1);
                d2= edst(X(K3(j),:),Mk2);
                d3= edst(X(K3(j),:),Mk3);
                d4= edst(X(K3(j),:),Mk4);
                dmin= min([d1,d2,d3,d4]);
                if (dmin == d1) && (length(K3)>1)
                    K1 = [K1 K3(j)];
                    K3 = [K3(1:j-1) K3(j+1:end)];
                    radi = true; j=j-1;
               
                elseif (dmin == d2) && (length(K3)>1)
                    K2 = [K2 K3(j)];
                    K3 = [K3(1:j-1) K3(j+1:end)];
                    radi = true; j=j-1;
                
                elseif (dmin == d4) && (length(K3)>1)
                    K4 = [K4 K3(j)];
                    K3 = [K3(1:j-1) K3(j+1:end)];
                    radi = true; j=j-1;
                end
                j=j+1;
            end
            
            j = 1;
            while (j <= length(K4))
                d1= edst(X(K4(j),:),Mk1);
                d2= edst(X(K4(j),:),Mk2);
                d3= edst(X(K4(j),:),Mk3);
                d4= edst(X(K4(j),:),Mk4);
                dmin= min([d1,d2,d3,d4]);
                if (dmin == d1) && (length(K4)>1)
                    K1 = [K1 K4(j)];
                    K4 = [K4(1:j-1) K4(j+1:end)];
                    radi = true; j=j-1;
                
                elseif (dmin == d2) && (length(K4)>1)
                    K2 = [K2 K4(j)];
                    K4 = [K4(1:j-1) K4(j+1:end)];
                    radi = true; j=j-1;
                
                elseif (dmin == d3) && (length(K4)>1)
                    K3 = [K3 K4(j)];
                    K4 = [K4(1:j-1) K4(j+1:end)];
                    radi = true; j=j-1;
                end
                j=j+1;
            end
            Mk1 = mean(X(K1,:),1);
            Mk2 = mean(X(K2,:),1);
            Mk3 = mean(X(K3,:),1);
            Mk4 = mean(X(K4,:),1);
        end
        itt_vect(d) = itt_vect(d)+i/Nitt;
    end
end


function v = ndek(X,M,S,P)
    v = 0.5*(X-M)'*(S^-1)*(X-M)+0.5*log(det(S))-0.5*log(P);
end

function d = edst(A,B)
    d = ((A(1)-B(1))^2+(A(2)-B(2))^2)^(1/2);
end