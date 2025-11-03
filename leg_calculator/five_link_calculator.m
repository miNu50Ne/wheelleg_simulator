clear;
tic;

syms phi1 phi4;
syms dphi1 dphi2 dphi3 dphi4;
syms l1 l2 l3 l4 l5;

x_B = l1 * cos(phi1);
y_B = l1 * sin(phi1);
x_D = l5 + l4 * cos(phi4);
y_D = l4 * sin(phi4);

l_BD = sqrt((x_D - x_B) * (x_D - x_B) + (y_D - y_B) * (y_D - y_B));
A0 = 2 * l2 * (x_D - x_B);
B0 = 2 * l2 * (y_D - y_B);
C0 = l2 * l2 + l_BD * l_BD - l3 * l3;

phi2 = 2 * atan2((B0 + sqrt(A0 * A0 + B0 * B0 - C0 * C0)), (A0 + C0));
phi3 = atan2(y_B - y_D + l2 * sin(phi2), x_B - x_D + l2 * cos(phi2));

x_C = l1 * cos(phi1) + l2 * cos(phi2);
y_C = l1 * sin(phi1) + l2 * sin(phi2);

L0 = sqrt((x_C - l5 / 2) * (x_C - l5 / 2) + y_C * y_C);
phi0 = atan2(y_C, x_C - l5 / 2);
theta = simplify(pi / 2 - phi0);

j11 = (l1 * sin(phi0 - phi3) * sin(phi1 - phi2)) / sin(phi3 - phi2);
j12 = (l4 * sin(phi0 - phi2) * sin(phi3 - phi4)) / sin(phi3 - phi2);
j21 = (l1 * cos(phi0 - phi3) * sin(phi1 - phi2)) / (L0 * sin(phi3 - phi2));
j22 = (l4 * cos(phi0 - phi2) * sin(phi3 - phi4)) / (L0 * sin(phi3 - phi2));
J = [j11 j12; j21 j22];

R = [cos(phi0 - pi / 2) -sin(phi0 - pi / 2);
     sin(phi0 - pi / 2) cos(phi0 - pi / 2)];
M = [0 -1 / L0;
     1 0];
M_v = [1 0;
       0 -1 / L0];

transform_matrix = J.' * R * M;
leg_posture = [L0; theta];

syms F Tp;
virtual_force = [F; Tp];
joint_torque = transform_matrix * virtual_force;

l1_ = 0.21; l2_ = 0.25; l5_ = 0; l3_ = l2_; l4_ = l1_; 
% [L0; theta] = leg_posture(phi1, phi4)
matlabFunction(leg_posture, 'File', '../function/leg_posture');
% [T1; T2] = joint_torque(F, Tp, phi1, phi4)
matlabFunction(joint_torque,'File','../function/joint_torque');
toc;
