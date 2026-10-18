%% Step 2c - Hill surfaces 2U(x,y,z) = C in 3D, Earth-Moon CR3BP
% Lower half (z <= 0) only, so the cut face shows the plane z = 0.
% References (course readers, Fundamentals of Astrodynamics, TU Delft):
%   K.F. Wakker, Fundamentals of Astrodynamics, TU Delft, 2015,
%     Sec. 3.6, Eq. (3.59).

clear; clc; close all
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')))

mu = 1/(1 + 81.30056);
[L, CL] = lagrangePoints(mu);
T  = plotTheme();

% 2U on a 3D grid (body at rest), lower half only

xv = linspace(-2, 2, 201);
yv = linspace(-2, 2, 201);
zv = linspace(-1.2, 0, 61);
[X, Y, Z] = meshgrid(xv, yv, zv);
P    = [X(:), Y(:), Z(:), zeros(numel(X), 3)];
TwoU = reshape(jacobiConstant(P, mu), size(X));
clear P                                       % ~120 MB, free it

C0     = jacobiConstant([0.8 0 0 0 0.3 0], mu);   % step 1 spacecraft
Cset   = [3.20, C0];
titles = {'C = 3.20 > C_{L1}: trapped', ...
          sprintf('C = %.3f (step 1): L1 and L2 open', C0)};

fig = figure('Color', T.bg, 'Position', [100 100 1300 600]);
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:2
    ax = nexttile(tl); hold(ax, 'on')
    C  = Cset(k);

    % Hill surface 2U = C
    
    s = patch(ax, isosurface(X, Y, Z, TwoU, C));
    isonormals(X, Y, Z, TwoU, s)
    set(s, 'FaceColor', T.accent1, 'EdgeColor', 'none', 'FaceAlpha', 0.6)

    % Fill the cut faces of the forbidden region (2U < C)
    
    c = patch(ax, isocaps(X, Y, Z, TwoU, C, 'below'));
    set(c, 'FaceColor', T.body, 'EdgeColor', 'none')

    % Primaries and libration points in the plane z = 0
    
    zm = 0.02;                                    % just above the lid
    plot3(ax, -mu,  0, zm, 'o', 'MarkerSize', 8, 'MarkerFaceColor', T.accent2, 'MarkerEdgeColor', 'none')
    plot3(ax, 1-mu, 0, zm, 'o', 'MarkerSize', 4, 'MarkerFaceColor', T.fg,      'MarkerEdgeColor', 'none')
    plot3(ax, L(:,1), L(:,2), zm*ones(5,1), '+', 'Color', T.accent3, 'MarkerSize', 7, 'LineWidth', 1.2)
    styleAxes(ax, T)
    set(ax, 'ZColor', T.fg)
       axis(ax, 'equal'); xlim(ax, [-2 2]); ylim(ax, [-2 2]); zlim(ax, [-1.2 0.05])
    view(ax, -35, 35)
    camlight(ax); lighting(ax, 'gouraud'); material(ax, 'dull')
    xlabel(ax, 'x'); ylabel(ax, 'y'); zlabel(ax, 'z')
    title(ax, titles{k})
end

title(tl, 'Hill surfaces 2U(x,y,z) = C, cut at z = 0', 'Color', T.fg, 'FontName', T.font)
exportgraphics(fig, fullfile(root, 'figures', 'step2_hillSurface.png'), ...
               'Resolution', 200, 'BackgroundColor', T.bg)