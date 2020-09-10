clear variables; close all; clc;

%% generisanje podataka
N = 500; % br odbiraka

P11 = 0.6;
P12 = 0.4;
M11 = [1; 1];
M12 = [6; 4];
S11 = [4, 1.1; 1.1, 2];
S12 = [3, -0.8; -0.8, 1.5];

P21 = 0.55;
P22 = 0.45;
M21 = [7; -4];
M22 = [6; 0];
S21 = [2, 1.1; 1.1, 4];
S22 = [3, 0.8; 0.8, 0.5];

X1 = zeros(1,N);
X2 = zeros(2,N);

    