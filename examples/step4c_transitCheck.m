%% Step 4c-3 - The tube explains the transit window of step 3
% Cuts the Earth-side stable tube with the section x = 0.80 (step 3 launch
% line), predicts the edges of the transit window from where the cut crosses
% y = 0, and checks them by brute force.
% Transit = reaching x = 1 - mu (Moon realm) before coming back to x = 0.80.
% References:
%   Koon, Lo, Marsden & Ross, Dynamical Systems, the Three-Body Problem
%     and Space Mission Design, 2nd ed., 2022, Sec. 2.7 (p. 52, tubes as
%     separatrices), Sec. 2.9 (Conley's wedge), Sec. 4.4 (Steps 4-5).

clear; clc; close all
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')))

load(fullfile(root, 'data', 'lyapunovL1_C318.mat'), 'orbit', 'mu', 'Ct')
L   = lagrangePoints(mu);
T   = plotTheme();
day = 27.321661/(2*pi);
xE  = 0.80;  xM = 1 - mu;
V0  = sqrt(jacobiConstant([xE 0 0 0 0 0], mu) - Ct);     % step 3 launch speed

% 1) Cut of the Earth-side stable tube with x = 0.80

s0 = [orbit(2); 0; 0; 0; orbit(3); 0];
N  = 400;
[~, UsE] = lyapunovManifoldICs(s0, orbit(4), mu, 'stable', N);
opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-12, 'Events', @(t,s) section(t, s, xE));
cut = nan(N, 6);

for k = 1:N

    [~, ~, ~, ye] = ode113(@(t,s) cr3bpEOM(t, s, mu), [0 -30], UsE(k,:).', opts);
    
    if ~isempty(ye), cut(k,:) = ye(1,:); end
end

fprintf('Tube cut: %d of %d trajectories reached x = %.2f, C on cut = %.10f\n', ...
        nnz(~isnan(cut(:,1))), N, xE, mean(jacobiConstant(cut, mu)));

% 2) Window edges: where the closed curve (y, vy) crosses y = 0

y  = cut(:,2);  vy = cut(:,5);
j  = [2:N, 1].';                               % next point (closed curve)
iz = find(y.*y(j) < 0);
vyc = vy(iz) + (vy(j(iz)) - vy(iz)).*(y(iz)./(y(iz) - y(j(iz))));
thEdge = sort(asind(vyc/V0));
fprintf('Tube prediction: transit for %.2f < theta < %.2f deg\n', thEdge(1), thEdge(end));

% 3) Brute force: launch from (0.80, 0) every 1 deg, as in step 3

th  = -90:1:90;
tM  = nan(size(th));
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @(t,s) transitEvent(t, s, xE, xM));

for k = 1:numel(th)

    sk = [xE; 0; 0; V0*cosd(th(k)); V0*sind(th(k)); 0];
    [~, ~, te, ~, ie] = ode113(@(t,s) cr3bpEOM(t, s, mu), [0 12], sk, opts);
    
    if ~isempty(ie) && ie(end) == 1, tM(k) = te(end); end
end

isT = ~isnan(tM);
fprintf('Brute force:     transit for %d <= theta <= %d deg (1 deg steps)\n', ...
        min(th(isT)), max(th(isT)));
th5 = -180:5:175;  in5 = th5 > thEdge(1) & th5 < thEdge(end);
fprintf('Step 3 grid (5 deg): %d directions inside the tube: %s deg\n', nnz(in5), mat2str(th5(in5)));

% ---------------------------------------------------------------------

fig = figure('Color', T.bg, 'Position', [100 100 1250 560]);
tl  = tiledlayout(fig, 1, 2, 'Padding', 'loose', 'TileSpacing', 'compact');

% (a) Poincare section x = 0.80: tube cut and step 3 launch states

ax1 = nexttile(tl); hold(ax1, 'on')
fill(ax1, y, vy, T.accent4, 'FaceAlpha', 0.15, 'EdgeColor', T.accent4, 'LineWidth', 1.5)
plot(ax1, zeros(1, nnz(~isT)), V0*sind(th(~isT)), '.', 'Color', T.grid, 'MarkerSize', 10)
plot(ax1, zeros(1, nnz(isT)),  V0*sind(th(isT)),  '.', 'Color', T.accent2, 'MarkerSize', 10)
styleAxes(ax1, T)
xlabel(ax1, 'y'); ylabel(ax1, 'v_y')
legend(ax1, {'Stable tube cut', 'Launch: no transit', 'Launch: transit'}, 'TextColor', T.fg, ...
       'Color', T.bg, 'EdgeColor', T.grid, 'Location', 'southwest')
title(ax1, sprintf('Section x = %.2f (step 3 launch line)', xE))

% (b) time to reach the Moon realm against launch direction

ax2 = nexttile(tl); hold(ax2, 'on')
plot(ax2, th(isT), tM(isT)*day, 'o-', 'Color', T.accent2, 'MarkerFaceColor', T.accent2, 'MarkerSize', 4)
xline(ax2, thEdge(1), '--', 'Color', T.accent4, 'LineWidth', 1.2)
xline(ax2, thEdge(end), '--', 'Color', T.accent4, 'LineWidth', 1.2)
styleAxes(ax2, T)
xlim(ax2, [-30 90])
xlabel(ax2, 'Launch direction \theta [deg]'); ylabel(ax2, 'Time to reach x = 1 - \mu [days]')
title(ax2, 'Brute force (dots) vs tube prediction (dashed)')

title(tl, 'Step 4c-3: the stable tube is the transit window', 'Color', T.fg, 'FontName', T.font)
set(fig, 'InvertHardcopy', 'off')
print(fig, fullfile(root, 'figures', 'step4c_transitCheck.png'), '-dpng', '-r200')

% ---------------------------------------------------------------------
function [val, isterminal, direction] = section(~, s, xs)

    val = s(1) - xs;  isterminal = 1;  direction = 0;

end

function [val, isterminal, direction] = transitEvent(~, s, xE, xM)

%TRANSITEVENT  1) reaches the Moon line x = xM (transit)
%              2) comes back to x = xE moving towards the Earth (no transit)

    val        = [s(1) - xM; s(1) - xE];
    isterminal = [1; 1];
    direction  = [1; -1];

end