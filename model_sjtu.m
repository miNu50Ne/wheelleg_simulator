%Replica the model of Shanghai Jiaotong University

clear;
tic;

syms theta_w_l theta_w_r
syms d_theta_w_l d_theta_w_r
syms dd_theta_w_l dd_theta_w_r

syms s theta_l_l theta_l_r theta_b phi
syms d_s d_theta_l_l d_theta_l_r d_theta_b d_phi
syms dd_s dd_theta_l_l dd_theta_l_r dd_theta_b dd_phi

syms T_w_l T_w_r T_b_l T_b_r

syms R_w R_l l_l l_r
syms l_w_l l_w_r
syms l_b_l l_b_r
syms l_c

syms m_w m_l m_b
syms I_w I_l_l I_l_r I_b I_z

%parameters
Rw = 0.07;
Rl = 0.22;

eta_l = 0.2945;
lwl = eta_l * l_l + 0.0368;
lbl = l_l - lwl;
lbr = lbl;
lwr = lwl;

lc = 0;

mw = 0.2;
mb = 16.98;
ml = 0.86;

Iw = 0.00036;
Ib = 0.36;
Iz = 0.33;
Ill = 0.3404 * l_l * l_l + 0.0463 * l_l + 0.0161;
Ilr = 0.3404 * l_r * l_r + 0.0463 * l_r + 0.0161;

g = 9.80665;

eq_1 = (I_w * l_l / R_w + m_w * R_w * l_w_l + m_l * R_w * l_b_l) * dd_theta_w_l ...
    + (m_l * l_w_l * l_b_l - I_l_l) * dd_theta_l_l ...
    + (m_l * l_w_l + 0.5 * m_b * l_b_l) * 9.81 * theta_l_l ...
    + T_b_l - T_w_l * (1 + l_l / R_w) == 0;

eq_2 = (I_w * l_r / R_w + m_w * R_w * l_w_r + m_l * R_w * l_b_r) * dd_theta_w_r ...
    + (m_l * l_w_r * l_b_r - I_l_r) * dd_theta_l_r ...
    + (m_l * l_w_r + 0.5 * m_b * l_b_r) * 9.81 * theta_l_r ...
    + T_b_r - T_w_r * (1 + l_r / R_w) == 0;

eq_3 =- (m_w * R_w ^ 2 + I_w + m_l * R_w ^ 2 + 0.5 * m_b * R_w ^ 2) * dd_theta_w_l ...
    - (m_w * R_w ^ 2 + I_w + m_l * R_w ^ 2 + 0.5 * m_b * R_w ^ 2) * dd_theta_w_r ...
    - (m_l * R_w * l_w_l + 0.5 * m_b * R_w * l_w_l) * dd_theta_l_l ...
    - (m_l * R_w * l_w_r + 0.5 * m_b * R_w * l_w_r) * dd_theta_l_r ...
    + T_w_l + T_w_r == 0;

eq_4 = (m_w * R_w * l_c + I_w * l_c / R_w + m_l * R_w * l_c) * dd_theta_w_l ...
    + (m_w * R_w * l_c + I_w * l_c / R_w + m_l * R_w * l_c) * dd_theta_w_r ...
    + m_l * l_w_l * l_c * dd_theta_l_l ...
    + m_l * l_w_r * l_c * dd_theta_l_r ...
    - I_b * dd_theta_b ...
    + m_b * 9.81 * l_c * theta_b ...
    - (T_w_l + T_w_r) * l_c / R_w ...
    - (T_b_l + T_b_r) == 0;

eq_5 = (0.5 * I_z * R_w / R_l + I_w * R_l / R_w) * dd_theta_w_l ...
    - (0.5 * I_z * R_w / R_l + I_w * R_l / R_w) * dd_theta_w_r ...
    + 0.5 * I_z * l_l / R_l * dd_theta_l_l ...
    - 0.5 * I_z * l_r / R_l * dd_theta_l_r ...
    - T_w_l * R_l / R_w + T_w_r * R_l / R_w == 0;

[dd_theta_l_l, dd_theta_l_r, dd_theta_w_l, dd_theta_w_r, dd_theta_b] = solve([eq_1, eq_2, eq_3, eq_4, eq_5], [dd_theta_l_l, dd_theta_l_r, dd_theta_w_l, dd_theta_w_r, dd_theta_b]);

dd_s = 0.5 * R_w * (dd_theta_w_l + dd_theta_w_r);

dd_phi = 0.5 * R_w / R_l * (-dd_theta_w_l + dd_theta_w_r) ...
    - 0.5 * l_l / R_l * cos(theta_l_l) * dd_theta_l_l ...
    + 0.5 * l_r / R_l * cos(theta_l_r) * dd_theta_l_r ...
    + 0.5 * l_l / R_l * sin(theta_l_l) * d_theta_l_l ^ 2 ...
    - 0.5 * l_r / R_w * sin(theta_l_r) * d_theta_l_r ^ 2;

from = [R_w, R_l, ...
            l_w_l, l_w_r, l_b_l, l_b_r, l_c, ...
            m_w, m_l, m_b, I_w, I_l_l, I_l_r, I_b, I_z];
to = [Rw, Rl, ...
          lwl, lwr, lbl, lbr, lc, ...
          mw, ml, mb, Iw, Ill, Ilr, Ib, Iz];

J_A = vpa(subs(subs(jacobian([d_s, dd_s, d_phi, dd_phi, d_theta_l_l, dd_theta_l_l, d_theta_l_r, dd_theta_l_r, d_theta_b, dd_theta_b], ...
    [s, d_s, phi, d_phi, theta_l_l, d_theta_l_l, theta_l_r, d_theta_l_r, theta_b, d_theta_b]), from, to), ...
    [s, d_s, phi, d_phi, theta_l_l, d_theta_l_l, theta_l_r, d_theta_l_r, theta_b, d_theta_b], ...
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
J_B = vpa(subs(subs(jacobian([d_s, dd_s, d_phi, dd_phi, d_theta_l_l, dd_theta_l_l, d_theta_l_r, dd_theta_l_r, d_theta_b, dd_theta_b], ...
    [T_w_l, T_w_r, T_b_l, T_b_r]), from, to), ...
    [theta_l_l, theta_l_r, T_w_l, T_w_r, T_b_l, T_b_r], ...
    [0, 0, 0, 0, 0, 0]));

%Functions of the leg length
A_calc = matlabFunction(J_A);
B_calc = matlabFunction(J_B);

%leg length
L_l_s = 0.1:0.01:0.45;
L_r_s = L_l_s;

K_s = zeros(4, 10, length(L_l_s), length(L_r_s));

% X = s d_s phi d_phi theta_l_l d_theta_l_l theta_l_r d_theta_l_r theta_b d_theta_b
Q = diag([500 100 50 5 1 1 1 1 1000 1]);
% U = T_w_l T_w_r T_b_l T_b_r
R = diag([1 1 0.25 0.25]);

T_s = 0.001;

for i = 1:length(L_l_s)

    for j = 1:length(L_r_s)
        A_c = A_calc(L_l_s(i), L_r_s(j));
        B_c = B_calc(L_l_s(i), L_r_s(j));

        [A_d, B_d] = c2d(A_c, B_c, T_s);

        K_s(:, :, i, j) = dlqr(A_d, B_d, Q, R);
    end

end

R_square = zeros(4, 10);
fit_type = fittype('poly22');
result = sym('result', [4, 10]);

K_result = sym('K', [4, 10]);
syms L;

[X_s, Y_s] = meshgrid(L_l_s, L_r_s);
x = X_s(:);
y = Y_s(:);

for i = 1:4

    for j = 1:10
        k = squeeze(K_s(i, j, :, :));
        k = k(:);

        [result, gof] = fit([x, y], k, fit_type);
        R_square(i, j) = gof.rsquare;

        K_result(i, j) = subs(subs(str2sym(formula(result)), ...
            coeffnames(result).', coeffvalues(result)), x, L);
    end

end

K_calc = matlabFunction(K_result);

toc;
