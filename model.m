clear;
tic;

% Displacement
syms theta(t) x(t) phi(t);
% Velocity
syms d_theta d_x d_phi;
% Acceleration
syms dd_theta dd_x dd_phi;
% Force
syms T T_p N_f N P N_M P_M;
% Inertia
syms m_w m_p M I_w I_p I_M;
% Parameters
syms l L L_M R

% Constants
syms g

N_M = M * diff((x + (L + L_M) * sin(theta) - l * sin(phi)), t, 2);
N = m_p * diff((x + L * sin(theta)), t, 2);
P_M = M * g + M * diff(((L + L_M) * cos(theta) + l * cos(phi)), t, 2);
P = P_M + m_p * g + m_p * diff((L * cos(theta)), t, 2);

eq_1 = diff(x, t, 2) == (T - N * R) / (I_w / R +m_w * R);
eq_2 = diff(theta, t, 2) == ((P * L + P_M * L_M) * sin(theta) - (N * L + N_M * L_M) * cos(theta) -T + T_p) / I_p;
eq_3 = diff(phi, t, 2) == (T_p + N_M * l * cos(phi) + P_M * l * sin(phi)) / I_M;

diffs = [diff(x(t), t, t), diff(theta(t), t, t), diff(phi(t), t, t), ...
             diff(x(t), t), diff(theta(t), t), diff(phi(t), t)];
vars = [dd_x, dd_theta, dd_phi, d_x, d_theta, d_phi];

eq_10 = subs(eq_1, diffs, vars);
eq_20 = subs(eq_2, diffs, vars);
eq_30 = subs(eq_3, diffs, vars);

[dd_x, dd_theta, dd_phi] = solve([eq_10, eq_20, eq_30], [dd_x, dd_theta, dd_phi]);

A = jacobian([d_theta, dd_theta, d_x, dd_x, d_phi, dd_phi], [theta, d_theta, x, d_x, phi, d_phi]);
B = jacobian([d_theta, dd_theta, d_x, dd_x, d_phi, dd_phi], [T, T_p]);

toc;
