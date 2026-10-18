%% Step 2d - Animation: zero-velocity curves as C decreases (necks opening)
% References (course readers, Fundamentals of Astrodynamics, TU Delft):
%   K.F. Wakker, Fundamentals of Astrodynamics, TU Delft, 2015,
%     Sec. 3.6, Eq. (3.59); Sec. 3.7, Eqs. (3.68)-(3.69).

clear; clc; close all
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')))

mu = 1/(1 + 81.30056);
[L, CL] = lagrangePoints(mu);
T  = plotTheme();

% 2U on the plane z = 0 (computed once, reused for every frame)

xv = linspace(-1.5, 1.5, 601);
[X, Y] = meshgrid(xv, xv);
TwoU = reshape(jacobiConstant([X(:), Y(:), zeros(numel(X), 4)], mu), size(X));

Cv      = linspace(3.22, 2.97, 120);              % decreasing C = more energy
gifFile = fullfile(root, 'figures', 'step2_necks.gif');

if isfile(gifFile), delete(gifFile), end

fig = figure('Color', T.bg, 'Position', [100 100 640 640]);
tl  = tiledlayout(fig, 1, 1, 'Padding', 'loose');   % margin so the title isn't clipped
ax  = nexttile(tl);

for k = 1:numel(Cv)
    C = Cv(k);
    cla(ax); hold(ax, 'on')
    zvcPanel(ax, T, xv, xv, X, Y, TwoU, C, mu, L)

    if     C > CL(1), state = 'trapped around Earth or Moon';
    elseif C > CL(2), state = 'L1 open: Earth-Moon transfer';
    elseif C > CL(3), state = 'L2 open: escape possible';
    elseif C > CL(4), state = 'L3 open';
    else,             state = 'no forbidden region in the plane';
    end

    title(ax, sprintf('C = %.4f   %s', C, state))
    drawnow

    % Capture the whole figure (no cropping) and append it to the GIF
    
    [A, map] = rgb2ind(frame2im(getframe(fig)), 256);
    dt = 0.08 + 1.5*(k == numel(Cv));             % pause 1.5 s on the last frame
    if k == 1
        imwrite(A, map, gifFile, 'gif', 'LoopCount', Inf, 'DelayTime', dt)
    else
        imwrite(A, map, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', dt)
    end

end