function P = Obelezja(X) 
%vraca jedan parametar P koji je dvodimenzionalan, a prima sliku

y=double(X); %posto imread vraca uint8, prebacujemo u double zbog racunskih operacija
% y=255-y; %negativ slike (ne mora da se radi)

%Binarizacija slike
T=0.8; %prag za binarizaciju  (0.2 ako je negativ)
y(y<T*max(max(y)))=0;
y(y>=T*max(max(y)))=255;

X=uint8(y); %vracamo u uint8 (255-y ako je negativ)
Z=X;



%Predobrada slike


[x,y]=size(X);
X=X(round(0.1*x):round(0.9*x),round(0.1*y):round(0.9*y));
[x,y]=size(X);
poc=1;
while (poc<x) && ((sum(X(poc,1:y))/y>252)...
   || (sum(X(poc,1:y))/y<=252 && sum(X(poc+5,1:y))/y>252))
    poc=poc+1;
end

kraj=x;
while (kraj>poc)&&((sum(X(kraj,1:y))/y>252)....
   || (sum(X(kraj,1:y))/y<=252 && sum(X(kraj-5,1:y))/y>252))
    kraj=kraj-1;
end

levo=1;
while (levo<y)&&((sum(X(1:x,levo))/x>252)...
  || (sum(X(1:x,levo))/x<=252 && sum(X(1:x,levo+5))/x>252))
    levo=levo+1;
end

desno=y;
while (desno>1)&&((sum(X(1:x,desno))/x>252)...
|| ((sum(X(1:x,desno))/x<=252) && sum(X(1:x,desno-5))/x>252))
    desno=desno-1;
end

X=X(poc:kraj,levo:desno); % izdvojena cifra
[x,y]=size(X);
% figure(1);
% subplot(2,1,1);
% imshow(X);
% subplot(2,1,2);
% imshow(Z);

%obelezje 1: posmatramo gustinu crnih piksela u pravougaoniku u sredini slike
P(1,1) = mean(X(round(2/5*x):round(3/5*x),round(1/3*y):round(2/3*y)), 'All');

cnt1 = 0;
for i=1:y
    ones = find(X(:,i)==0);
    if (ones)
        cnt1 = cnt1 + 1;
        for j=2:length(ones)
            if(ones(j)-ones(j-1)>10)
                cnt1 = cnt1 + 1;
            end
        end
    end
end
cnt2 = 0;
for i=1:round(x/2)
    ones = find(X(i,:)==0);
    if (ones)
        cnt2 = cnt2 + 1;
        for j=2:length(ones)
            if(ones(j)-ones(j-1)>10)
                cnt2 = cnt2 + 1;
            end
        end
    end
end
cnt3 = 0;
for i=round(x/2):x
    ones = find(X(i,:)==0);
    if (ones)
        cnt3 = cnt3 + 1;
        for j=2:length(ones)
            if(ones(j)-ones(j-1)>10)
                cnt3 = cnt3 + 1;
            end
        end
    end
end
%obelezje 2: odnos broja horizontalnih i vertikalnih linija
P(2,1) = (cnt1/y) / ((cnt2/x) + (cnt3/x));
%obelezje 3: broj vertikalnih linija u donjoj polovini
P(3,1) = (cnt3/x);
%obelezje 4: racuna gustinu piksela u pravougaoniku u donjem desnom uglu
%P(4,1) = mean(X(round(3/4*x):end,round(3/4*y):end),'All'); 
%P(4,1) = mean(mean(X(1:round(1/2*y),:))) - mean(mean(X(round(1/2*y):end)));
P(4,1) = mean(X(round(4/5*x):end,round(1/3*y):round(2/3*y)), 'All');
