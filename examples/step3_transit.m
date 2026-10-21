%% Step 3 - Transit through the L1 neck: which launch directions reach the Moon?
% Same C, same start point, same speed; only the launch direction changes.
% C_L2 < C < C_L1: only the L1 neck is open (Earth <-> Moon, no escape).
% References (course readers, Fundamentals of Astrodynamics, TU Delft):
%   K.F. Wakker, Fundamentals of Astrodynamics, TU Delft, 2015,
%     Sec. 3.3, Eq. (3.50); Sec. 3.4, Eq. (3.53); Secs. 3.6-3.7.

clear; clc; close all
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')))

mu  = 1/(1 + 81.30056);
[L, CL] = lagrangePoints(mu);
xL1 = L(1,1);
T   = plotTheme();

% Units: length = Earth-Moon distance, time = 1/omega

day   = 27.321661/(2*pi);                     % 1 time unit [days] (sidereal month/2pi)
vUnit = 384400/(day*86400);                   % 1 velocity unit [km/s]

% Initial conditions

C     = 3.18;                                 % C_L2 < C < C_L1
x0    = 0.80;                                 % Earth side of the neck
V0    = sqrt(jacobiConstant([x0 0 0 0 0 0], mu) - C);   % V^2 = 2U - C, Eq. (3.53)
theta = -180:5:175;                           % launch direction [deg] from +x
tEnd  = 6;                                    % ~26 days

% Integrate one trajectory per direction; event = first crossing of x = xL1

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,s) crossL1(t, s, xL1));
n      = numel(theta);
traj   = cell(n, 1);
tCross = nan(n, 1);

for k = 1:n
    th = deg2rad(theta(k));
    s0 = [x0; 0; 0; V0*cos(th); V0*sin(th); 0];
    [~, S, te] = ode113(@(t,s) cr3bpEOM(t, s, mu), [0 tEnd], s0, opts);
    traj{k} = S(:, 1:2);
    if ~isempty(te), tCross(k) = te(1); end
end

isT = ~isnan(tCross);

fprintf('C = %.4f  (C_L2 = %.4f, C_L1 = %.4f)\n', C, CL(2), CL(1));
fprintf('V0 = %.4f  = %.0f m/s in the rotating frame\n', V0, V0*vUnit*1e3);
fprintf('Transits in %.0f days: %d of %d directions\n', tEnd*day, nnz(isT), n);
fprintf('Window: %g to %g deg\n', min(theta(isT)), max(theta(isT)));
fprintf('Time to cross L1: %.1f to %.1f days\n', min(tCross)*day, max(tCross)*day);

% ---------------------------------------------------------------------

fig = figure('Color', T.bg, 'Position', [100 100 1250 560]);
tl  = tiledlayout(fig, 1, 2, 'Padding', 'loose', 'TileSpacing', 'compact');

% (a) trajectories near the neck

ax1 = nexttile(tl); hold(ax1, 'on')
xv = linspace(0.55, 1.25, 701);
yv = linspace(-0.35, 0.35, 351);
[X, Y] = meshgrid(xv, yv);
TwoU = reshape(jacobiConstant([X(:), Y(:), zeros(numel(X), 4)], mu), size(X));
zvcPanel(ax1, T, xv, yv, X, Y, TwoU, C, mu, L)
for k = find(~isT).'
    plot(ax1, traj{k}(:,1), traj{k}(:,2), 'Color', T.grid, 'LineWidth', 0.5)
end
for k = find(isT).'
    plot(ax1, traj{k}(:,1), traj{k}(:,2), 'Color', T.accent2, 'LineWidth', 0.9)
end
plot(ax1, L(1:2,1), L(1:2,2), '+', 'Color', T.accent3, 'MarkerSize', 8, 'LineWidth', 1.4)
plot(ax1, x0, 0, 'o', 'MarkerSize', 6, 'MarkerFaceColor', T.fg, 'MarkerEdgeColor', 'none')
xlim(ax1, [0.55 1.25]); ylim(ax1, [-0.35 0.35])
title(ax1, sprintf('C = %.2f: transit (orange) vs no transit (grey)', C))

% (b) transit window

ax2 = nexttile(tl); hold(ax2, 'on')
plot(ax2, theta(isT), tCross(isT)*day, 'o-', 'Color', T.accent2, ...
     'MarkerFaceColor', T.accent2, 'MarkerSize', 5)
plot(ax2, theta(~isT), zeros(nnz(~isT), 1), 'x', 'Color', T.grid, 'MarkerSize', 6)
styleAxes(ax2, T)
xlim(ax2, [-180 180]); xticks(ax2, -180:60:180)
xlabel(ax2, 'Launch direction \theta [deg]'); ylabel(ax2, 'Time to cross L1 [days]')
title(ax2, 'Transit window')

title(tl, 'Step 3: transit through the L1 neck, same energy, different direction', ...
      'Color', T.fg, 'FontName', T.font)
set(fig, 'InvertHardcopy', 'off')
print(fig, fullfile(root, 'figures', 'step3_transit.png'), '-dpng', '-r200')

% ---------------------------------------------------------------------

function [val, isterminal, direction] = crossL1(~, s, xL1)

%CROSSL1  Event: x crosses xL1 moving towards the Moon.
    val        = s(1) - xL1;    % zero when the spacecraft is at x = xL1
    isterminal = 0;             % keep integrating (we want the whole path)
    direction  = 1;             % only count crossings with x increasing
end