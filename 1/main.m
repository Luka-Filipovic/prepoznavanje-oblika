clear variables; close all; clc;

%% ucitavanje slika
baza = dir('BAZA/*.bmp');

N = numel(baza);
N1 = N/5;
NO = round(3/4*N1);
NT = N1 - NO;
Na = 0;
Ne = 0;
Ni = 0;
No = 0;
Nu = 0;
for i=1:N
    switch baza(i).name(5)
    case 'A'
        Na = Na + 1;
        Xa(Na,:) = Obelezja(imread(strcat("BAZA/", baza(i).name)));
    case 'E'
        Ne = Ne + 1;
	    Xe(Ne,:) = Obelezja(imread(strcat("BAZA/", baza(i).name)));
    case 'I'
        Ni = Ni + 1;
	    Xi(Ni,:) = Obelezja(imread(strcat("BAZA/", baza(i).name)));
    case 'O'
        No = No + 1;
        Xo(No,:) = Obelezja(imread(strcat("BAZA/", baza(i).name)));
    case 'U'
        Nu = Nu + 1;
        Xu(Nu,:) = Obelezja(imread(strcat("BAZA/", baza(i).name)));
    end
end
clearvars baza;
No = length(Xa(1,:)); %broj obelezja
%% ponavljanje cele klasifikacije Nitt puta
Nitt = 1000;
MK = zeros(5,5,Nitt); %vektor matrica konfuzije
for itt = 1:Nitt

%% nasumicno deljenje na testirajuci i obucavajuci skup
Xa=Xa(randperm(length(Xa)),:);
Xe=Xe(randperm(length(Xe)),:);
Xi=Xi(randperm(length(Xi)),:);
Xo=Xo(randperm(length(Xo)),:);
Xu=Xu(randperm(length(Xu)),:);

Xa_o=Xa(1:NO,:);
Xe_o=Xe(1:NO,:);
Xi_o=Xi(1:NO,:);
Xo_o=Xo(1:NO,:);
Xu_o=Xu(1:NO,:);

Xa_t=Xa(NO+1:end,:);
Xe_t=Xe(NO+1:end,:);
Xi_t=Xi(NO+1:end,:);
Xo_t=Xo(NO+1:end,:);
Xu_t=Xu(NO+1:end,:);

%% izracunavanje verovatnoce nad izabranim obucavajucim skupom
Sa = cov(Xa_o);
Se = cov(Xe_o);
Si = cov(Xi_o);
So = cov(Xo_o);
Su = cov(Xu_o);

Ma = mean(Xa_o)';
Me = mean(Xe_t)';
Mi = mean(Xi_t)';
Mo = mean(Xo_t)';
Mu = mean(Xu_t)';

%% klasifikacija i popunjavanje konfuzione matrice 
Mk = zeros(5,5); %konfuziona matrica
for k = 1:5
    
    switch(k) %odabiramo klasu za testiranje
        case 1
            Xt = Xa_t;
        case 2
            Xt = Xe_t;
        case 3
            Xt = Xi_t;
        case 4
            Xt = Xo_t;
        case 5
            Xt = Xu_t;
    end

    for i = 1:length(Xt)
        xt = Xt(i,:)';
        
        %racunamo fgv za svaku klasu sa izabranim parametrima
        fa = 1/(2*pi*det(Sa)^0.5)*exp(-0.5*(xt-Ma)'*Sa^(-1)*(xt-Ma));
        fe = 1/(2*pi*det(Se)^0.5)*exp(-0.5*(xt-Me)'*Se^(-1)*(xt-Me));
        fi = 1/(2*pi*det(Si)^0.5)*exp(-0.5*(xt-Mi)'*Si^(-1)*(xt-Mi));
        fo = 1/(2*pi*det(So)^0.5)*exp(-0.5*(xt-Mo)'*So^(-1)*(xt-Mo));
        fu = 1/(2*pi*det(Su)^0.5)*exp(-0.5*(xt-Mu)'*Su^(-1)*(xt-Mu));
        
        %odabiramo smestamo odbirak u klasu max fgv
        m = max([fa fe fi fo fu]);
        switch m
            case fa
                Mk(k,1) = Mk(k,1) + 1;
            case fe
                Mk(k,2) = Mk(k,2) + 1;
            case fi
                Mk(k,3) = Mk(k,3) + 1;
            case fo
                Mk(k,4) = Mk(k,4) + 1;
            case fu
                Mk(k,5) = Mk(k,5) + 1;
        end
    end
end
MK(:,:,itt) = Mk;
end

%% prikaz rezultata klasifikacije
disp("Primer matrice konfuzije:");
T = array2table(Mk,'VariableNames',{'A','E','I','O','U'},'RowName',...
    {'A  |','E  |','I  |','O  |','U  |'}); 
disp(T);
disp("Procenat tacnosti:");
disp(sum(diag(Mk))/(5*NT) * 100);
disp("Procenat greske:");
disp(100 - sum(diag(Mk))/(5*NT) * 100);

MK_mean = mean(MK,3);

disp("Usrednjena matrica konfuzije:");
T = array2table(MK_mean,'VariableNames',{'A','E','I','O','U'},'RowName',...
    {'A  |','E  |','I  |','O  |','U  |'}); 
disp(T);
disp("Usrednjen procenat tacnosti:");
disp(sum(diag(MK_mean))/(5*NT) * 100);
disp("Usrednjen procenat greske:");
disp(100-sum(diag(MK_mean))/(5*NT) * 100);
disp("Varijansa procenta greske:");
disp(var(100*(1-reshape([MK(1,1,:);MK(2,2,:);MK(3,3,:);MK(4,4,:);MK(5,5,:)],[],1)/(5*NT))));