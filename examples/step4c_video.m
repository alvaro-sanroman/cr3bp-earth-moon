%% STEP 4c-4 - Cinematic 3D video: transit through the L1 neck
% A spacecraft with Jacobi constant C = 3.18 starts in the Earth realm,
% rides the stable manifold tube of the L1 Lyapunov orbit through the
% neck and ends captured in the Moon realm (the L2 neck is closed).
%
% Scene (rotating frame, nondimensional units, z = 0 is the orbital plane):
%   - Earth and Moon as plain lit spheres, drawn at 2x their true radius
%   - forbidden region (2U < C) and zero-velocity curve on the z = 0 plane
%   - Lyapunov orbit around L1 (white) and its manifold tubes
%       stable   (green)  : branches integrated backward in time
%       unstable (orange) : branches integrated forward in time
%       Moon branches up to x = 1 - mu, Earth branches up to x = 0.80
%   - spacecraft: launched from (0.80, 0) at theta = 28 deg, |v| from C,
%     integrated from -8 to +25 days (t = 0 is the crossing of x = 0.80)
%
% Sources: Wakker, Fundamentals of Astrodynamics (TU Delft course reader),
% Ch. 3; Koon, Lo, Marsden & Ross, Dynamical Systems, the Three-Body
% Problem and Space Mission Design (2022), Sec. 4.4 (Eq. 4.4.2-4.4.3).
%
% Output: figures/cr3bp_transit3d.mp4  (+ a still: figures/cr3bp_transit3d_still.png)
% Keep the figure window visible while it renders (do not minimise,
% dock or close it): getframe captures what is on screen.

clear; close all; clc
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')))
T = plotTheme();

%% 1. Data from step 4c-1
S   = load(fullfile(root, 'data', 'lyapunovL1_C318.mat'));
mu  = S.mu;
x0  = S.orbit(2);
vy0 = S.orbit(3);
Tp  = S.orbit(4);
C   = S.orbit(5);

[L, ~] = lagrangePoints(mu);
if size(L, 1) < size(L, 2), L = L.'; end    % one row per point
xE  = -mu;            % Earth
xM  = 1 - mu;         % Moon
xL1 = L(1, 1);
xL2 = L(2, 1);

TU  = 27.321661/(2*pi);              % time unit [days]
RE  = 2*6378/384400;                 % Earth radius x2 (nondimensional)
RM  = 2*1737/384400;                 % Moon radius x2  (nondimensional)
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

%% 2. Lyapunov orbit and its monodromy matrix
w0 = [x0; 0; 0; 0; vy0; 0; reshape(eye(6), [], 1)];
[~, Wend] = ode113(@(t, w) cr3bpEOMstm(t, w, mu), [0 Tp], w0, opts);
if norm(Wend(end, 1:6)' - w0(1:6)) > 1e-6     % saved T was a half period
    Tp = 2*Tp;
end

N = 24;                                      % points on the orbit
tN = linspace(0, Tp, N + 1);
[~, W] = ode113(@(t, w) cr3bpEOMstm(t, w, mu), tN, w0, opts);
M = reshape(W(end, 7:42), 6, 6);             % monodromy matrix

[V, D] = eig(M);
lam = diag(D);
[~, iu] = max(abs(lam));
[~, is] = min(abs(lam));
Yu0 = real(V(:, iu));  Yu0 = Yu0/norm(Yu0);  if Yu0(1) < 0, Yu0 = -Yu0; end
Ys0 = real(V(:, is));  Ys0 = Ys0/norm(Ys0);  if Ys0(1) < 0, Ys0 = -Ys0; end

[~, Sorb] = ode113(@(t, s) cr3bpEOM(t, s, mu), linspace(0, Tp, 400), w0(1:6), opts);

%% 3. Manifold tubes (KLMR Eq. 4.4.2-4.4.3)
epsil = 1e-6;
tubes = struct('xy', {}, 'type', {});        % type 1 = unstable, 2 = stable
for k = 1:N
    s   = W(k, 1:6)';
    Phi = reshape(W(k, 7:42), 6, 6);
    Yu  = Phi*Yu0;  Yu = Yu/norm(Yu);
    Ys  = Phi*Ys0;  Ys = Ys/norm(Ys);

    % {initial state, time direction, stop at x =, type}
    cases = { s + epsil*Yu,  +10, xM,   1;    % unstable, Moon side
              s - epsil*Yu,  +10, 0.80, 1;    % unstable, Earth side
              s + epsil*Ys,  -10, xM,   2;    % stable,   Moon side
              s - epsil*Ys,  -10, 0.80, 2 };  % stable,   Earth side
    for c = 1:4
        oe = odeset(opts, 'Events', @(t, q) xSection(t, q, cases{c, 3}));
        [~, Q] = ode113(@(t, q) cr3bpEOM(t, q, mu), [0 cases{c, 2}], cases{c, 1}, oe);
        tubes(end + 1).xy = Q(:, 1:2);       %#ok<SAGROW>
        tubes(end).type   = cases{c, 4};
    end
end

%% 4. Spacecraft trajectory, theta = 28 deg from (0.80, 0)
TwoU = @(x, y) x.^2 + y.^2 + 2*(1 - mu)./hypot(x + mu, y) + 2*mu./hypot(x - 1 + mu, y);
V0   = sqrt(TwoU(0.80, 0) - C);
th   = 28;
s0   = [0.80; 0; 0; V0*cosd(th); V0*sind(th); 0];
solF = ode113(@(t, s) cr3bpEOM(t, s, mu), [0  25/TU], s0, opts);
solB = ode113(@(t, s) cr3bpEOM(t, s, mu), [0  -8/TU], s0, opts);

tGrid = linspace(-8, 25, 8000);              % days
P = zeros(3, numel(tGrid));
iB = tGrid < 0;
P(:, iB)  = deval(solB, tGrid(iB)/TU, 1:3);
P(:, ~iB) = deval(solF, tGrid(~iB)/TU, 1:3);

tCrossL1 = tGrid(find(P(1, :) > xL1, 1));    % ~ 1.54 days
tCrossM  = tGrid(find(P(1, :) > xM,  1));    % ~ 4.78 days

%% 5. Figure and scene
fig = figure('Color', T.bg, 'Units', 'pixels', 'Position', [60 60 1600 900], ...
    'WindowStyle', 'normal', 'Resize', 'off', 'MenuBar', 'none', ...
    'ToolBar', 'none', 'InvertHardcopy', 'off', 'Name', 'CR3BP transit 3D');
ax = axes(fig, 'Units', 'normalized', 'Position', [0 0 1 1], 'Color', T.bg);
hold(ax, 'on')
axis(ax, 'off')
set(ax, 'DataAspectRatio', [1 1 1], 'Clipping', 'off', 'SortMethod', 'depth', ...
    'XLim', [-9 9], 'YLim', [-9 9], 'ZLim', [-9 9], 'Projection', 'perspective')

blend = @(col, a) T.bg + a*(col - T.bg);     % faded colour without transparency

% Stars (far away, for parallax)
rng(7)
nS = 700;
dS = randn(3, nS);  dS = dS./vecnorm(dS);
dS = 8*dS;
scatter3(ax, dS(1, :), dS(2, :), dS(3, :), 3 + 6*rand(1, nS), ...
    'filled', 'MarkerFaceColor', T.fg, 'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor', 'none');

% Forbidden region and zero-velocity curve on the plane z = 0
xv = linspace(-1.5, 1.5, 500);
[Xg, Yg] = meshgrid(xv, xv);
TU2 = TwoU(Xg, Yg);
forbid = double(TU2 < C);
surf(ax, Xg, Yg, -0.002*ones(size(Xg)), 'FaceColor', blend(T.fg, 0.18), ...
    'EdgeColor', 'none', 'FaceAlpha', 'interp', 'AlphaData', 0.55*forbid, ...
    'AlphaDataMapping', 'none', 'FaceLighting', 'none');
Cc = contourc(xv, xv, TU2, [C C]);
k = 1;
while k < size(Cc, 2)
    n = Cc(2, k);
    plot3(ax, Cc(1, k+1:k+n), Cc(2, k+1:k+n), zeros(1, n), ...
        'Color', blend(T.fg, 0.45), 'LineWidth', 0.8);
    k = k + n + 1;
end

% Faint reference grid: circles around the barycentre and the Earth-Moon line
for r = [0.5 1.0 1.5]
    a = linspace(0, 2*pi, 300);
    plot3(ax, r*cos(a), r*sin(a), -0.001*ones(size(a)), 'Color', blend(T.fg, 0.12), 'LineWidth', 0.6);
end
plot3(ax, [-1.5 1.5], [0 0], [-0.001 -0.001], '--', 'Color', blend(T.fg, 0.15), 'LineWidth', 0.6);

% Earth and Moon: plain lit spheres
[sx, sy, sz] = sphere(80);
surf(ax, xE + RE*sx, RE*sy, RE*sz, 'FaceColor', [0.30 0.50 0.85], ...
    'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'AmbientStrength', 0.35, ...
    'DiffuseStrength', 0.9, 'SpecularStrength', 0.2);
surf(ax, xM + RM*sx, RM*sy, RM*sz, 'FaceColor', [0.70 0.70 0.72], ...
    'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'AmbientStrength', 0.35, ...
    'DiffuseStrength', 0.9, 'SpecularStrength', 0.05);
light('Parent', ax, 'Style', 'infinite', 'Position', [-1 -0.7 0.9], 'Color', [1 1 1]);
light('Parent', ax, 'Style', 'infinite', 'Position', [1 0.5 0.4], 'Color', 0.25*[1 1 1]);

text(ax, xE, 0, RE + 0.05, 'Earth', 'Color', T.fg, 'FontName', T.font, ...
    'FontSize', 13, 'HorizontalAlignment', 'center', 'Clipping', 'off');
text(ax, xM, 0, RM + 0.03, 'Moon', 'Color', T.fg, 'FontName', T.font, ...
    'FontSize', 13, 'HorizontalAlignment', 'center', 'Clipping', 'off');

% L1 and L2
plot3(ax, [xL1 xL2], [0 0], [0.003 0.003], 'd', 'MarkerSize', 6, ...
    'MarkerFaceColor', T.fg, 'MarkerEdgeColor', 'none');
text(ax, xL1, 0, 0.035, 'L_1', 'Color', T.fg, 'FontName', T.font, 'FontSize', 13, ...
    'HorizontalAlignment', 'center', 'Clipping', 'off');
text(ax, xL2, 0, 0.035, 'L_2', 'Color', T.fg, 'FontName', T.font, 'FontSize', 13, ...
    'HorizontalAlignment', 'center', 'Clipping', 'off');

% Manifold tubes (start invisible, fade in during the intro)
colU = T.accent2;                            % unstable: orange
colS = T.accent4;                            % stable:   green
hTube = gobjects(numel(tubes), 1);
for k = 1:numel(tubes)
    hTube(k) = plot3(ax, tubes(k).xy(:, 1), tubes(k).xy(:, 2), zeros(size(tubes(k).xy, 1), 1), ...
        'Color', T.bg, 'LineWidth', 0.9);
end
tubeCol = {colU, colS};
tubeAlpha = 0.55;

% Lyapunov orbit
plot3(ax, Sorb(:, 1), Sorb(:, 2), 0.001*ones(size(Sorb, 1), 1), 'Color', T.fg, 'LineWidth', 2.2);

% Spacecraft: trail + glowing head
colSC = T.accent1;
hTrail = plot3(ax, nan, nan, nan, 'Color', colSC, 'LineWidth', 2.0);
p0 = scPos(-8, solB, solF, TU);
hGlow2 = scatter3(ax, p0(1), p0(2), p0(3), 520, 'filled', 'MarkerFaceColor', colSC, ...
    'MarkerFaceAlpha', 0.10, 'MarkerEdgeColor', 'none');
hGlow1 = scatter3(ax, p0(1), p0(2), p0(3), 170, 'filled', 'MarkerFaceColor', colSC, ...
    'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
hHead  = scatter3(ax, p0(1), p0(2), p0(3), 36, 'filled', 'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', 'none');

%% 6. Captions (figure coordinates, independent of the camera)
txt = @(pos, str, fs, varargin) annotation(fig, 'textbox', pos, 'String', str, ...
    'Color', T.fg, 'EdgeColor', 'none', 'FontName', T.font, 'FontSize', fs, ...
    'Interpreter', 'tex', 'FitBoxToText', 'off', 'VerticalAlignment', 'top', varargin{:});

txt([0.03 0.955 0.70 0.04], '\bfEarth-Moon CR3BP: transit through the L_1 neck', 24);
txt([0.03 0.905 0.80 0.04], sprintf(['Rotating frame   |   \\mu = %.5f   |   C = %.2f   |   ' ...
    'Lyapunov orbit around L_1 and its invariant manifold tubes'], mu, C), 13, ...
    'Color', blend(T.fg, 0.75));
hTime = txt([0.70 0.955 0.27 0.04], '', 18, 'HorizontalAlignment', 'right');

rgb = @(c) sprintf('%.3f,%.3f,%.3f', c);
keyStr = { ['\color[rgb]{' rgb(T.fg)  '}\bf---\rm  Lyapunov orbit (L_1)'], ...
           ['\color[rgb]{' rgb(colS)  '}\bf---\rm  stable manifold tube'], ...
           ['\color[rgb]{' rgb(colU)  '}\bf---\rm  unstable manifold tube'], ...
           ['\color[rgb]{' rgb(colSC) '}\bf---\rm  spacecraft, C = 3.18'], ...
           ['\color[rgb]{' rgb(blend(T.fg, 0.6)) '}\bf---\rm  forbidden region (2U < C)'] };
txt([0.78 0.235 0.21 0.18], keyStr, 12);
txt([0.03 0.045 0.94 0.04], sprintf(['Earth and Moon drawn at 2x true radius.  Spacecraft launched ' ...
    'from (x, y) = (0.80, 0) at \\theta = %d\\circ with |v| fixed by C;  t = 0 at that point.  ' ...
    'Tubes: N = %d trajectories per branch, \\epsilon = 10^{-6}.'], th, N), 10.5, ...
    'Color', blend(T.fg, 0.6));

%% 7. Camera keyframes (u = progress of the main sequence, 0..1)
%        u     target            dist  az    el
key = [ 0.00   0.30   0  0      3.0   -110  35
        0.30   xL1    0  0      0.6    -70  30
        0.45   0.90   0  0      0.5    -40  28
        0.75   xM     0  0      0.7      0  35
        1.00   0.40   0  0      2.6     40  40 ];
camVA = 38;                                       % view angle [deg]
camAt = @(u) deal(  [pchip(key(:,1), key(:,2), u), pchip(key(:,1), key(:,3), u), pchip(key(:,1), key(:,4), u)], ...
                    exp(pchip(key(:,1), log(key(:,5)), u)), ...
                    pchip(key(:,1), key(:,6), u), ...
                    pchip(key(:,1), key(:,7), u));
% Time warp: more frames near the neck
tOfU = @(u) pchip([0 0.2 0.6 1], [-8 0 5 25], u);

%% 8. Video
outFile = fullfile(root, 'figures', 'cr3bp_transit3d.mp4');
vw = VideoWriter(outFile, 'MPEG-4');
vw.FrameRate = 30;
vw.Quality   = 95;
open(vw);
cleanupVideo = onCleanup(@() close(vw));

nIntro = 75;  nMain = 750;  nOutro = 75;
setCam(ax, [0.3 0 0], 3.0, -150, 35, camVA);
drawnow
f0 = getframe(fig);
H = 2*floor(size(f0.cdata, 1)/2);
Wd = 2*floor(size(f0.cdata, 2)/2);

% Intro: tubes fade in while the camera swings round (no captions)
for f = 1:nIntro
    s = f/nIntro;  e = 0.5 - 0.5*cos(pi*s);       % ease in/out
    setCam(ax, [0.3 0 0], 3.0, -150 + 40*e, 35, camVA);
    for k = 1:numel(tubes)
        hTube(k).Color = blend(tubeCol{tubes(k).type}, tubeAlpha*e);
    end
    hTime.String = sprintf('t = %+5.1f days', -8);
    drawnow
    writeVideo(vw, grabFrame(fig, H, Wd));
end

% Main: spacecraft flies, camera follows the keyframes
stillDone = false;
for f = 1:nMain
    u  = (f - 1)/(nMain - 1);
    td = tOfU(u);
    [tg, d, az, el] = camAt(u);
    setCam(ax, tg, d, az, el, camVA);

    p = scPos(td, solB, solF, TU);
    in = tGrid <= td;
    set(hTrail, 'XData', [P(1, in) p(1)], 'YData', [P(2, in) p(2)], 'ZData', [P(3, in) p(3)] + 0.002);
    set([hHead hGlow1 hGlow2], 'XData', p(1), 'YData', p(2), 'ZData', p(3) + 0.002);

    hTime.String = sprintf('t = %+5.1f days', td);
    drawnow
    fr = grabFrame(fig, H, Wd);
    writeVideo(vw, fr);

    if ~stillDone && td > 2.2                        % still for the README
        imwrite(fr.cdata, fullfile(root, 'figures', 'cr3bp_transit3d_still.png'));
        stillDone = true;
    end
end

% Outro: slow pull-back with the full trajectory on screen
[tg, d, az, el] = camAt(1);
for f = 1:nOutro
    s = f/nOutro;  e = 0.5 - 0.5*cos(pi*s);
    setCam(ax, tg, d*(1 + 0.1*e), az + 20*e, el, camVA);
    drawnow
    writeVideo(vw, grabFrame(fig, H, Wd));
end

clear cleanupVideo                                    % closes the MP4
fprintf('Video saved: %s  (%d frames, %.1f s)\n', outFile, ...
    nIntro + nMain + nOutro, (nIntro + nMain + nOutro)/30);
fprintf('Crosses x_L1 at %.3f days, x = 1-mu at %.3f days\n', tCrossL1, tCrossM);

%% Local functions
function p = scPos(td, solB, solF, TU)
    if td < 0
        p = deval(solB, td/TU, 1:3);
    else
        p = deval(solF, td/TU, 1:3);
    end
end

function setCam(ax, tgt, d, az, el, va)
    dirv = [cosd(el)*cosd(az), cosd(el)*sind(az), sind(el)];
    set(ax, 'CameraTarget', tgt, 'CameraPosition', tgt + d*dirv, ...
        'CameraUpVector', [0 0 1], 'CameraViewAngle', va);
end

function fr = grabFrame(fig, H, W)
    fr = getframe(fig);
    c = fr.cdata;
    if size(c, 1) ~= H || size(c, 2) ~= W
        c2 = zeros(H, W, 3, 'uint8');
        h = min(H, size(c, 1));  w = min(W, size(c, 2));
        c2(1:h, 1:w, :) = c(1:h, 1:w, :);
        c = c2;
    end
    fr.cdata = c;
    fr.colormap = [];
end

function [val, isterminal, direction] = xSection(~, s, xs)
    val = s(1) - xs;
    isterminal = 1;
    direction = 0;
end