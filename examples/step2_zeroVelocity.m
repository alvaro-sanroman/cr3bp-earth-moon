%% Step 2b - Zero-velocity curves (Hill regions), Earth-Moon CR3BP
% Forbidden region: 2U < C, since V^2 = 2U - C >= 0.
% References (course readers, Fundamentals of Astrodynamics, TU Delft):
%   K.F. Wakker, Fundamentals of Astrodynamics, TU Delft, 2015,
%     Sec. 3.4, Eq. (3.53); Sec. 3.6, Eq. (3.59); Sec. 3.7.
clear; clc; close all
root = fileparts(fileparts(mfilename('fullpath')));   % project folder
addpath(genpath(fullfile(root, 'src')))

mu = 1/(1 + 81.30056);
[L, CL] = lagrangePoints(mu);

% 2U on a grid in the plane z = 0 (= C of a body at rest)
xv = linspace(-1.5, 1.5, 801);
yv = linspace(-1.5, 1.5, 801);
[X, Y] = meshgrid(xv, yv);
P    = [X(:), Y(:), zeros(numel(X), 4)];
TwoU = reshape(jacobiConstant(P, mu), size(X));

% Step 1 trajectory (same initial state)
s0   = [0.8; 0; 0; 0; 0.3; 0];
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[~, S] = ode113(@(t,s) cr3bpEOM(t,s,mu), [0 30], s0, opts);
C0 = jacobiConstant(s0.', mu);
V  = sqrt(sum(S(:,4:6).^2, 2));

fprintf('C0 = %.6f  (C_L1 = %.6f, C_L2 = %.6f, C_L3 = %.6f, C_L4 = %.6f)\n', C0, CL(1:4));
fprintf('Min speed along trajectory: %.2e  (cusps on the curve)\n', min(V));

% One panel per regime: C between consecutive C_L values
Cset   = [C0, 3.20, 3.18, 3.10, 3.00, 2.98];
titles = {sprintf('Step 1 spacecraft, C = %.3f', C0), ...
          'C > C_{L1}', 'C_{L2} < C < C_{L1}', 'C_{L3} < C < C_{L2}', ...
          'C_{L4} < C < C_{L3}', 'C < C_{L4}'};

T   = plotTheme();
fig = figure('Color', T.bg, 'Position', [100 100 1200 760]);
tl  = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:6
    ax = nexttile(tl); hold(ax, 'on')
    zvcPanel(ax, T, xv, yv, X, Y, TwoU, Cset(k), mu, L)
    if k == 1
        plot(ax, S(:,1), S(:,2), 'Color', T.accent2, 'LineWidth', 0.8)
    end
    title(ax, titles{k})
end
title(tl, 'Zero-velocity curves 2U = C, Earth-Moon system', ...
      'Color', T.fg, 'FontName', T.font)
exportgraphics(fig, fullfile(root, 'figures', 'step2_zeroVelocity.png'), ...
               'Resolution', 200, 'BackgroundColor', T.bg)

