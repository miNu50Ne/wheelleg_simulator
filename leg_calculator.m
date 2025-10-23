clear;
tic;
%The geometry part
syms phi_1 phi_2;
syms l_1 l_2 l_3 l_4 l_5;
syms L_0 phi_0;

phi_3 = 0.5 * (pi - phi_1 -phi_2);
phi_4 = acos(l_1 / l_2 * cos(phi_3));
alpha = pi - 2 * phi_4;
beta = pi - phi_2 - alpha;
phi_5 = pi - beta;

x_B = l_1 * cos(phi_1);
x_G = (l_1 + l_4) * cos(phi_2);
y_G = (l_1 + l_4) * sin(phi_2);
x_H = x_G + l_5 * cos(phi_5);
y_H = y_G + l_5 * sin(phi_5);
