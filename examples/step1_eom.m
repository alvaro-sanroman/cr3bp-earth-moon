%% STEP 1 - CR3BP equations of motion and Jacobi constant check
%  Integrates one trajectory in the Earth-Moon rotating frame and checks
%  that Jacobi's constant stays constant (accuracy of the integration).
%  Ref: Wakker, Fundamentals of Astrodynamics (TU Delft), Eqs. (3.50), (3.53).

clear; clc; close all
root   = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
figDir = fullfile(root, 'figures');

% -------------------------------------------------------- parameters
% Earth/Moon mass ratio 81.30056 (IAU 2009) -> mu = m2/(m1+m2)

mu = 1/(1 + 81.30056);

% initial state [x y z vx vy vz] (non-dimensional, rotating frame)

s0    = [0.8; 0; 0; 0; 0.3; 0];
tSpan = [0 30];               % 30 time units ~ 130 days (1 unit ~ 4.35 d)

% ------------------------------------------------------- integration

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t, S] = ode113(@(t, s) cr3bpEOM(t, s, mu), tSpan, s0, opts);

% ------------------------------------------------ Jacobi constant check

C     = jacobiConstant(S, mu);
drift = abs(C - C(1));

fprintf('mu                 = %.9f\n', mu);
fprintf('Jacobi constant C0 = %.6f\n', C(1));
fprintf('max |C - C0|       = %.2e\n', max(drift));
fprintf('steps taken        = %d\n', numel(t));

% ------------------------------------------------------------- plot

T   = plotTheme();
fig = figure('Color', T.bg, 'Position', [100 100 1200 500], 'WindowStyle', 'normal');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) trajectory in the rotating frame

ax1 = nexttile(tl); hold(ax1, 'on');
plot(ax1, S(:,1), S(:,2), '-', 'Color', T.accent1, 'LineWidth', 0.8);
plot(ax1, -mu,   0, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0.30 0.55 1.00], 'MarkerEdgeColor', 'none');
plot(ax1, 1-mu,  0, 'o', 'MarkerSize', 6,  'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', 'none');
plot(ax1, s0(1), s0(2), 'o', 'MarkerSize', 6, 'MarkerFaceColor', T.accent2, 'MarkerEdgeColor', 'none');
axis(ax1, 'equal');
xlabel(ax1, 'x'); ylabel(ax1, 'y');
title(ax1, sprintf('Trajectory in the rotating frame, C = %.4f', C(1)));
lg = legend(ax1, {'trajectory', 'Earth', 'Moon', 'start'}, 'Location', 'northeast');
set(lg, 'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid);
styleAxes(ax1, T);

% (b) Jacobi constant error

ax2 = nexttile(tl);
semilogy(ax2, t, drift + eps, '-', 'Color', T.accent3, 'LineWidth', 1.2);
xlabel(ax2, 't (non-dimensional)'); ylabel(ax2, '|C - C_0|');
title(ax2, 'Jacobi constant: integration error');
styleAxes(ax2, T);

ax1.Toolbar.Visible = 'off';  ax2.Toolbar.Visible = 'off';
exportgraphics(fig, fullfile(figDir, 'step1_eom.png'), ...
               'Resolution', 200, 'BackgroundColor', 'current');