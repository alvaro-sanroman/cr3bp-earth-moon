%% Step 4c-2 - Invariant manifolds (tubes) of the L1 Lyapunov orbit at C = 3.18
% Stable and unstable manifolds, Earth and Moon branches, integrated until
% the sections x = 0.80 (Earth side, step 3 launch line) and x = 1 - mu (Moon).
% References:
%   Koon, Lo, Marsden & Ross, Dynamical Systems, the Three-Body Problem
%     and Space Mission Design, 2nd ed., 2022, Sec. 4.4, Step 4 (pp. 114-116).
%   K.F. Wakker, Fundamentals of Astrodynamics, TU Delft, 2015, Sec. 3.12.

clear; clc; close all
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')))

load(fullfile(root, 'data', 'lyapunovL1_C318.mat'), 'orbit', 'mu', 'Ct')
L   = lagrangePoints(mu);
T   = plotTheme();
day = 27.321661/(2*pi);

s0 = [orbit(2); 0; 0; 0; orbit(3); 0];        % Lyapunov orbit at C = Ct
Tp = orbit(4);
N  = 40;                                      % trajectories per branch
xE = 0.80;                                    % Earth-side section (step 3)
xM = 1 - mu;                                  % Moon-side section (KLMR U3 line)

% Initial states on both manifolds

[UuM, UuE] = lyapunovManifoldICs(s0, Tp, mu, 'unstable', N);
[UsM, UsE] = lyapunovManifoldICs(s0, Tp, mu, 'stable',   N);

% Branch: initial states, time direction, stop section, colour

br = {UuM, +1, xM, T.accent3;   UuE, +1, xE, T.accent3;
      UsM, -1, xM, T.accent4;   UsE, -1, xE, T.accent4};
names = {'Unstable, Moon side', 'Unstable, Earth side', 'Stable, Moon side', 'Stable, Earth side'};

fig = figure('Color', T.bg, 'Position', [100 100 1150 720]);
tl  = tiledlayout(fig, 1, 1, 'Padding', 'loose');
ax  = nexttile(tl); hold(ax, 'on')
xv = linspace(0.76, 1.02, 701);  yv = linspace(-0.11, 0.11, 593);
[X, Y] = meshgrid(xv, yv);
TwoU = reshape(jacobiConstant([X(:), Y(:), zeros(numel(X), 4)], mu), size(X));
zvcPanel(ax, T, xv, yv, X, Y, TwoU, Ct, mu, L)

h = gobjects(4, 1);

for b = 1:4

    X0 = br{b,1};  dir = br{b,2};  xs = br{b,3};
    opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-12, 'Events', @(t,s) section(t, s, xs));
    tHit = nan(N, 1);
  
    for k = 1:N

        [~, S, te] = ode113(@(t,s) cr3bpEOM(t, s, mu), [0 dir*30], X0(k,:).', opts);
        hk = plot(ax, S(:,1), S(:,2), 'Color', br{b,4}, 'LineWidth', 0.5);
       
        if ~isempty(te), tHit(k) = abs(te(1)); end

    end

    h(b) = hk;
    fprintf('%-22s reached the section: %2d of %d, time %.1f-%.1f days\n', ...
            names{b}, nnz(~isnan(tHit)), N, min(tHit)*day, max(tHit)*day);

end

% Lyapunov orbit, sections and step 3 launch point

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[~, S] = ode113(@(t,s) cr3bpEOM(t, s, mu), [0 Tp], s0, opts);
hL = plot(ax, S(:,1), S(:,2), 'Color', T.fg, 'LineWidth', 2);
xline(ax, xE, '--', 'Color', T.grid);  xline(ax, xM, '--', 'Color', T.grid);
plot(ax, xE, 0, 'o', 'MarkerSize', 6, 'MarkerFaceColor', T.accent2, 'MarkerEdgeColor', 'none')

xlim(ax, [0.76 1.02]); ylim(ax, [-0.11 0.11])
legend(ax, [hL h(1) h(3)], {'Lyapunov orbit (C = 3.18)', 'Unstable manifold', 'Stable manifold'}, ...
       'TextColor', T.fg, 'Color', T.bg, 'EdgeColor', T.grid, 'Location', 'southeast')
title(ax, 'Invariant manifold tubes of the L1 Lyapunov orbit')
title(tl, 'Step 4c-2: stable and unstable manifolds at C = 3.18', 'Color', T.fg, 'FontName', T.font)
set(fig, 'InvertHardcopy', 'off')
print(fig, fullfile(root, 'figures', 'step4c_manifolds.png'), '-dpng', '-r200')

% ---------------------------------------------------------------------

function [val, isterminal, direction] = section(~, s, xs)

%SECTION  Event: the trajectory reaches the line x = xs (either direction).

    val = s(1) - xs;  isterminal = 1;  direction = 0;
end