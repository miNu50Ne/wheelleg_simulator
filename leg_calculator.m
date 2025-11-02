clear;
tic;

%The geometry part
syms phi_1(t) phi_2(t);
syms l_1 l_2 l_3 l_4 l_5;
syms d_phi_1 d_phi_2;

syms F T_p;
syms phi_1_ phi_2_;

assumeAlso(l_1 < l_2);

phi_3 = 0.5 * (pi - phi_1 -phi_2);
phi_4 = acos(l_1 / l_2 * cos(phi_3));
alpha = pi - 2 * phi_4;
beta = pi - phi_2 - alpha;
phi_5 = pi - beta;

x_G = (l_1 + l_4) * cos(phi_2);
y_G = (l_1 + l_4) * sin(phi_2);
x_H = x_G + l_5 * cos(phi_5);
y_H = y_G + l_5 * sin(phi_5);

L_0 = sqrt(x_H ^ 2 + y_H ^ 2);
phi = atan2(y_H, x_H);
theta = simplify(pi/2 - phi);

x_dot_H = diff(x_H, t);
y_dot_H = diff(y_H, t);

x_dot_H = subs(x_dot_H, [diff(phi_1, t), diff(phi_2, t)], [d_phi_1, d_phi_2]);
y_dot_H = subs(y_dot_H, [diff(phi_1, t), diff(phi_2, t)], [d_phi_1, d_phi_2]);

x_dot = [x_dot_H; y_dot_H];
q_dot = [d_phi_1; d_phi_2];

x_dot = simplify(collect(x_dot, q_dot));

J = simplify(jacobian(x_dot, q_dot));

R = [sin(phi), cos(phi);
     -cos(phi), sin(phi)];

M = [0, -1 / L_0;
     1, 0];
torque_transform_matrix = simplify(J.' * R * M);

virtual_force = [F; T_p];
torque_transform_matrix = formula(subs(torque_transform_matrix, [phi_1, phi_2], [phi_1_, phi_2_]));
torque = torque_transform_matrix * virtual_force;

attitude = formula(subs([L_0, theta], [phi_1, phi_2], [phi_1_, phi_2_]));

% velocity;

% user parameter
l1_ = 0.0945; l2_ = 0.1125; l3_ = 0.065; l4_ = 0.1155; l5_ = 0.25;

torque = subs(torque, [l_1, l_2, l_3, l_4, l_5], [l1_, l2_, l3_, l4_, l5_]);
attitude = subs(attitude, [l_1, l_2, l_3, l_4, l_5], [l1_, l2_, l3_, l4_, l5_]);
matlabFunction(torque, 'File', 'torque_transform');
matlabFunction(attitude, 'File', 'leg_attitude');
toc;
