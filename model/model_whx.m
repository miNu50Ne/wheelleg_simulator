%Replica the model of Wang Hongxi

clear;
tic;

% Displacement
syms theta(t) x(t) phi(t);
syms theta0 x0 phi0;
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
syms leg_len
% Constants
syms g

R_1 = 0.09;
eta_l = 0.3088;
L_1 = eta_l * leg_len + 0.0380;
L_M_1 = leg_len - L_1;
l_1 = 0.035;
m_w_1 = 0.39;
m_p_1 = 1.26;
M_1 = 13;
I_w_1 = 0.005;
I_p_1 = 0.3202 * leg_len * leg_len + 0.0556 * leg_len + 0.0240;
I_M_1 = 0.2106;
g_0 = 9.80665;

N_M = M * diff((x + (L + L_M) * sin(theta) - l * sin(phi)), t, 2);
N = m_p * diff((x + L * sin(theta)), t, 2) + N_M;
P_M = M * g + M * diff(((L + L_M) * cos(theta) + l * cos(phi)), t, 2);
P = P_M + m_p * g + m_p * diff((L * cos(theta)), t, 2);

eq_1 = diff(x, t, 2) == (T - N * R) / (I_w / R +m_w * R);
eq_2 = diff(theta, t, 2) == ((P * L + P_M * L_M) * sin(theta) - (N * L + N_M * L_M) * cos(theta) -T + T_p) / I_p;
eq_3 = diff(phi, t, 2) == (T_p + N_M * l * cos(phi) + P_M * l * sin(phi)) / I_M;

from = [diff(x(t), t, t), diff(theta(t), t, t), diff(phi(t), t, t), ...
            diff(x(t), t), diff(theta(t), t), diff(phi(t), t), x, theta, phi];
to = [dd_x, dd_theta, dd_phi, d_x, d_theta, d_phi, x0, theta0, phi0];

eq_10 = subs(eq_1, from, to);
eq_20 = subs(eq_2, from, to);
eq_30 = subs(eq_3, from, to);

[dd_x, dd_theta, dd_phi] = solve([eq_10, eq_20, eq_30], [dd_x, dd_theta, dd_phi]);

A = subs(jacobian([d_theta, dd_theta, d_x, dd_x, d_phi, dd_phi], [theta0, d_theta, x0, d_x, phi0, d_phi]), [theta0, d_theta, d_x, phi0, d_phi, T, T_p], [0, 0, 0, 0, 0, 0, 0]);
A = subs(A, [R, L, L_M, l, m_w, m_p, M, I_w, I_p, I_M, g], [R_1, L_1, L_M_1, l_1, m_w_1, m_p_1, M_1, I_w_1, I_p_1, I_M_1, g_0]);
A = vpa(A);
B = subs(jacobian([d_theta, dd_theta, d_x, dd_x, d_phi, dd_phi], [T, T_p]), [theta0, d_theta, d_x, phi0, d_phi, T, T_p], [0, 0, 0, 0, 0, 0, 0]);
B = subs(B, [R, L, L_M, l, m_w, m_p, M, I_w, I_p, I_M, g], [R_1, L_1, L_M_1, l_1, m_w_1, m_p_1, M_1, I_w_1, I_p_1, I_M_1, g_0]);
B = vpa(B);

%Functions of the leg length
A_calc = matlabFunction(A);
B_calc = matlabFunction(B);

%leg length
L_0 = 0.1:0.01:0.45;
K = zeros(2, 6, length(L_0));

% X = theta d_theta x d_x phi d_phi
Q = diag([1, 1, 500, 100, 5000, 1]);
% U = T T_p
R = diag([1, 0.25]);

T_s = 0.001;

for index = 1:length(L_0)
    A_0 = A_calc(L_0(index));
    B_0 = B_calc(L_0(index));

    [A_d, B_d] = c2d(A_0, B_0, T_s);

    K(:, :, index) = dlqr(A_d, B_d, Q, R);
    % K(:, :, index) = lqr(A_0, B_0, Q, R);
end

R_square = zeros(2, 6);
fit_type = fittype('poly4');
result = sym('result', [2, 6]);

K_result = sym('K', [2, 6]);
syms L;

figure('Name', 'Fitted Curve', 'NumberTitle', 'off');
tiledlayout(2, 6, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:2

    for j = 1:6
        k = squeeze(K(i, j, :));
        [result, gof] = fit(L_0(:), k, fit_type);
        % Coefficient of determination -> goodness of fit
        R_square(i, j) = gof.rsquare;

        K_result(i, j) = subs(subs(str2sym(formula(result)), ...
            coeffnames(result).', coeffvalues(result)), x, L);

        % 绘图区域
        nexttile((i - 1) * 6 + j);
        hold on;

        % 绘制原始点
        plot(L_0, k, 'bo', 'MarkerSize', 2.5, 'DisplayName', '数据');

        % 拟合曲线
        fitted_vals = result(L_0);
        plot(L_0, fitted_vals, 'r-', 'LineWidth', 1.2, 'DisplayName', '拟合');

        % 标签和标题
        xlabel('L_0');
        ylabel(sprintf('K(%d,%d)', i, j));
        title(sprintf('R^2 = %.4f', R_square(i, j)));

        legend('Location', 'best');
        grid on;
        hold off;

    end

end

K_calc = matlabFunction(K_result);

L0 = 0.2;
K_0 = K_calc(L0);

disp("L=0.2:");
disp(K_0);

toc;
