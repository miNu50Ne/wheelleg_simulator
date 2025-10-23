clear;
tic;
%The geometry part
syms phi_1(t) phi_2(t);
syms l_1 l_2 l_3 l_4 l_5;
syms d_phi_1 d_phi_2;
syms L_0 phi_0;

phi_3 = 0.5 * (pi - phi_1 -phi_2);
phi_4 = acos(l_1 / l_2 * cos(phi_3));
alpha = pi - 2 * phi_4;
beta = pi - phi_2 - alpha;
phi_5 = pi - beta;

x_G = (l_1 + l_4) * cos(phi_2);
y_G = (l_1 + l_4) * sin(phi_2);
x_H = x_G + l_5 * cos(phi_5);
y_H = y_G + l_5 * sin(phi_5);

L_0 = sqrt(x_H^2 + y_H^2);
phi_0 = atan2(y_H,x_H);

x_dot_H = diff(x_H, t);
y_dot_H = diff(y_H, t);

x_dot_H = subs(x_dot_H,[diff(phi_1,t),diff(phi_2,t)],[d_phi_1,d_phi_2]);
y_dot_H = subs(x_dot_H,[diff(phi_1,t),diff(phi_2,t)],[d_phi_1,d_phi_2]);

x_dot = [x_dot_H; y_dot_H];
q_dot = [d_phi_1; d_phi_2];

x_dot = simplify(collect(x_dot,q_dot));
J = simplify(jacobian(x_dot,q_dot));

R = [sin(phi_0),cos(phi_0); -cos(phi_0),sin(phi_0)];

M = [0,-1/L_0;1,0];
torque_transform_matrix=simplify(J.'*R*M);
