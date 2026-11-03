%% Step 4b - Planar Lyapunov orbit around L1 by differential correction
% Seed from the linear solution, Newton correction of vy0 using the STM,
% then checks: period vs linear theory, Jacobi constant, periodicity,
% monodromy eigenvalues.
% References:
%   K.F. Wakker, Fundamentals of Astrodynamics, TU Delft, 2015,
%     Sec. 3.10, Eqs. (3.114)-(3.119) (course reader).
%   Koon, Lo, Marsden & Ross, Dynamical Systems, the Three-Body Problem
%     and Space Mission Design, 2nd ed., 2022, Secs. 2.7, 4.2-4.4,
%     Eqs. (2.7.2), (4.2.5), (4.4.1).

clear; clc; close all
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')))

mu  = 1/(1 + 81.30056);
L   = lagrangePoints(mu);
xe  = L(1,1);
T   = plotTheme();
day = 27.321661/(2*pi);                       % 1 time unit [days]
km  = 384400;                                 % 1 length unit [km]

% Linear constants at L1 (KLMR Eq. 2.6.2, p. 47, Eq. 2.7.2)

mb  = mu/abs(xe - 1 + mu)^3 + (1 - mu)/abs(xe + mu)^3;
lam = sqrt((mb - 2 + sqrt(9*mb^2 - 8*mb))/2);
nu  = sqrt(-(mb - 2 - sqrt(9*mb^2 - 8*mb))/2);
tau = -(nu^2 + 2*mb + 1)/(2*nu);

% Seed from the linear solution, KLMR Eq. (4.4.1)

Ax   = 1e-3;                                  % x-amplitude (~384 km)
x0   = xe - Ax;
vy0g = -Ax*nu*tau;

% Differential correction

tol = 1e-10;
[vy0, Tp, hist] = lyapunovCorrect(x0, vy0g, mu, tol);

% One full period with the STM -> monodromy matrix M = Phi(T)

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
w0   = [x0; 0; 0; 0; vy0; 0; reshape(eye(6), 36, 1)];
[~, W] = ode113(@(t,w) cr3bpEOMstm(t, w, mu), [0 Tp], w0, opts);
M    = reshape(W(end, 7:42), 6, 6);
lm   = eig(M);
[~, iS] = sort(abs(lm));  lm = lm(iS);
C    = jacobiConstant(W(:, 1:6), mu);

fprintf('Seed (Eq. 4.4.1):  x0 = %.8f   vy0 = %.8f\n', x0, vy0g);
fprintf('Corrected:                       vy0 = %.8f\n', vy0);
fprintf('|vx| at half period: '); fprintf('%.1e  ', hist); fprintf('\n');
fprintf('Period: %.6f = %.3f days   (linear 2pi/nu = %.3f days)\n', Tp, Tp*day, 2*pi/nu*day);
fprintf('C = %.6f   (C_L1 = %.6f)   max |C - C0| = %.1e\n', C(1), jacobiConstant([xe 0 0 0 0 0], mu), max(abs(C - C(1))));
fprintf('Periodicity error |s(T) - s(0)| = %.1e\n', max(abs(W(end,1:6) - W(1,1:6))));
fprintf('Monodromy eigenvalues:\n'); disp(lm)
fprintf('lambda1*lambda2 = %.6f   exp(lambda*T) = %.1f\n', abs(lm(1)*lm(end)), exp(lam*Tp));

% ---------------------------------------------------------------------

fig = figure('Color', T.bg, 'Position', [100 100 1250 560]);
tl  = tiledlayout(fig, 1, 2, 'Padding', 'loose', 'TileSpacing', 'compact');

% (a) corrected orbit vs linear ellipse

ax1 = nexttile(tl); hold(ax1, 'on')
tt  = linspace(0, 2*pi/nu, 300);
plot(ax1, -Ax*cos(nu*tt)*km, -Ax*tau*sin(nu*tt)*km, '--', 'Color', T.fg, 'LineWidth', 1)
plot(ax1, (W(:,1) - xe)*km, W(:,2)*km, 'Color', T.accent2, 'LineWidth', 2)
plot(ax1, 0, 0, '+', 'Color', T.accent3, 'MarkerSize', 9, 'LineWidth', 1.4)
styleAxes(ax1, T); axis(ax1, 'equal')
xlabel(ax1, 'x - x_{L1} [km]'); ylabel(ax1, 'y [km]')
legend(ax1, {'Linear (Eq. 4.4.1)', 'Corrected (nonlinear)', 'L1'}, ...
       'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid, 'Location', 'northeast')
title(ax1, sprintf('Lyapunov orbit, A_x = %.0f km, T = %.2f days', Ax*km, Tp*day))

% (b) Newton convergence

ax2 = nexttile(tl); hold(ax2, 'on')
semilogy(ax2, 1:numel(hist), hist, 'o-', 'Color', T.accent2, 'MarkerFaceColor', T.accent2)
yline(ax2, tol, '--', 'tolerance', 'Color', T.grid, 'LabelHorizontalAlignment', 'left')
styleAxes(ax2, T); set(ax2, 'YScale', 'log')
xlabel(ax2, 'Iteration'); ylabel(ax2, '|v_x| at the x-axis crossing')
title(ax2, 'Differential correction (Newton)')

title(tl, 'Step 4b: Lyapunov orbit around L1', 'Color', T.fg, 'FontName', T.font)
set(fig, 'InvertHardcopy', 'off')
print(fig, fullfile(root, 'figures', 'step4b_lyapunov.png'), '-dpng', '-r200')