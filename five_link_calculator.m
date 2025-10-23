% VMC solution

clear;
tic;
%The geometry part
syms phi_1(t) phi_2(t) phi_3(t) phi_4(t);
syms d_phi_1 d_phi_4;
syms l_1 l_2 l_3 l_4 l_5;
syms L_0 phi_0;

x_B = l_1 * cos(phi_1);
y_B = l_1 * sin(phi_1);
x_D = l_5 + l_4 * cos(phi_4);
y_D = l_4 * sin(phi_4);

x_C = x_B + l_2 * cos(phi_2);
y_C = y_B + l_2 * sin(phi_2);

x_dot_B = diff(x_B, t);
y_dot_B = diff(y_B, t);
x_dot_C = diff(x_C, t);
y_dot_C = diff(y_C, t);
x_dot_D = diff(x_D, t);
y_dot_D = diff(y_D, t);

d_phi_2 = ((x_dot_D - x_dot_B) * cos(phi_3) + (y_dot_D - y_dot_B) * sin(phi_3)) / (l_2 * sin(phi_3 - phi_2));

x_dot_C = subs(x_dot_C, diff(phi_2, t), d_phi_2);
x_dot_C = subs(x_dot_C, ...
    [diff(phi_1, t), diff(phi_4, t)], ...
    [d_phi_1, d_phi_4]);
y_dot_C = subs(y_dot_C, diff(phi_2, t), d_phi_2);
y_dot_C = subs(y_dot_C, ...
    [diff(phi_1, t), diff(phi_4, t)], ...
    [d_phi_1, d_phi_4]);

x_dot = [x_dot_C; y_dot_C];
q_dot = [d_phi_1; d_phi_4];
x_dot = simplify(collect(x_dot, q_dot));
%Jacobian matrix
J = simplify(jacobian(x_dot, q_dot));

%Rotation
R = [cos(phi_0), sin(phi_0);
     -sin(phi_0), cos(phi_0)];
%Transformation
M = [0, -1 / L_0;
     1, 0];
% %
% M_v = [1, 0;
%        0, -1 / L_0];

%Dynamics Jacobian matrix
torque_transform_matrix = simplify(J.' * R * M);

% %Inverse dynamics Jacobian matrix
% syms T_1 T_2

% %force_transform_matrix
% force_transform_matrix = simplify(torque_transform_matrix.' \ [T_1; T_2]);

% A_0 = 2 * l_2 * (x_D - x_B);
% B_0 = 2 * l_2 * (y_D - y_B);
% C_0 = l_2 ^ 2 +sqrt((x_D - x_B) ^ 2 + (y_D - y_B) ^ 2) - l_3 ^ 2;

% phi_2_ = simplify(2 * atan2(B_0 + sqrt(A_0 ^ 2 + B_0 ^ 2 - C_0 ^ 2), A_0 + C_0));
% phi_3_ = simplify(atan2(x_B - x_D + l_2 * cos(phi_2_), y_B - y_D + l_2 * sin(phi_2_)));

% x_C_ = x_B + l_2 * cos(phi_2_);
% y_C_ = y_B + l_2 * sin(phi_2_);

% L_0_ = simplify(sqrt((x_C - l_5 / 2) ^ 2 + y_C ^ 2));
% phi_0_ = atan2(y_C, x_C - l_5 / 2);

% syms phi_1_ phi_4_;
% syms d_phi_1_ d_phi_4_;
% syms F T_p;

% %leg posture
% leg_posture = subs(formula([L_0_, phi_0_]), [phi_1, phi_4], [phi_1_, phi_4_]);

% torque_transform_matrix = simplify(subs(torque_transform_matrix, [phi_2(t), phi_3(t), L_0, phi_0], [phi_2_, phi_3_, L_0_, phi_0_]));
% %Joint motor torque
% joint_torque = torque_transform_matrix * [F; T_p];
toc;
