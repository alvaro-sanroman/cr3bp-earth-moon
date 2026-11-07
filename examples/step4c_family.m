%% Step 4c-1 - Lyapunov family around L1 by natural-parameter continuation
% Builds the planar Lyapunov family from small amplitudes up to the orbit
% with C = 3.18 (the energy of step 3), and saves that orbit for step 4c-2.
% References:
%   Koon, Lo, Marsden & Ross, Dynamical Systems, the Three-Body Problem
%     and Space Mission Design, 2nd ed., 2022, Sec. 4.4, Steps 2-3
%     (Eq. 4.4.1, continuation p. 114), Sec. 4.3 (monodromy matrix).

clear; clc; close all
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')))

mu  = 1/(1 + 81.30056);
[L, CL] = lagrangePoints(mu);
xe  = L(1,1);
T   = plotTheme();
day = 27.321661/(2*pi);                       % 1 time unit [days]
km  = 384400;                                 % 1 length unit [km]

% Linear constants at L1 (KLMR Eq. 2.6.2, p. 47, Eq. 2.7.2)

mb  = mu/abs(xe - 1 + mu)^3 + (1 - mu)/abs(xe + mu)^3;
nu  = sqrt(-(mb - 2 - sqrt(9*mb^2 - 8*mb))/2);
tau = -(nu^2 + 2*mb + 1)/(2*nu);

Ct  = 3.18;                                   % target Jacobi constant (step 3)
dAx = 1e-3;                                   % continuation step (~384 km)

% ---------------------------------------------------------------------
% 1) Two small orbits from the linear seed, KLMR Eq. (4.4.1)
%    fam columns: [Ax, x0, vy0, T, C]

fam = zeros(2, 5);
for k = 1:2

    Ax = k*dAx;
    x0 = xe - Ax;
    [vy0, Tk] = lyapunovCorrect(x0, -Ax*nu*tau, mu);
    fam(k,:) = [Ax, x0, vy0, Tk, jacobiConstant([x0 0 0 0 vy0 0], mu)];

end

% 2) Continuation: extrapolate vy0 from the last two orbits (KLMR p. 114)

while fam(end,5) > Ct

    Ax    = fam(end,1) + dAx;
    slope = (fam(end,3) - fam(end-1,3)) / (fam(end,1) - fam(end-1,1));
    vyg   = fam(end,3) + slope*dAx;
    x0    = xe - Ax;
    [vy0, Tk] = lyapunovCorrect(x0, vyg, mu);
    fam(end+1,:) = [Ax, x0, vy0, Tk, jacobiConstant([x0 0 0 0 vy0 0], mu)]; %#ok<AGROW> 

end

% 3) Refine to C = Ct: secant method on Ax between the last two orbits

a = fam(end-1,:);  b = fam(end,:);

for k = 1:10

    Ax  = a(1) + (Ct - a(5))*(b(1) - a(1))/(b(5) - a(5));
    vyg = a(3) + (b(3) - a(3))*(Ax - a(1))/(b(1) - a(1));
    x0  = xe - Ax;
    [vy0, Tk] = lyapunovCorrect(x0, vyg, mu);
    a = b;
    b = [Ax, x0, vy0, Tk, jacobiConstant([x0 0 0 0 vy0 0], mu)];
   
    if abs(b(5) - Ct) < 1e-10, break, end

end

orbit = b;                                    % [Ax, x0, vy0, T, C] at C = Ct

% 4) Monodromy matrix of the target orbit

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
s0   = [orbit(2); 0; 0; 0; orbit(3); 0];
[~, W] = ode113(@(t,w) cr3bpEOMstm(t, w, mu), [0 orbit(4)], [s0; reshape(eye(6), 36, 1)], opts);
M    = reshape(W(end, 7:42), 6, 6);
lm   = eig(M);  [~, iS] = sort(abs(lm));  lm = lm(iS);

% ---------------------------------------------------------------------

fprintf('   Ax [km]     T [days]     C\n');
fprintf('   %6.0f      %7.3f     %.6f\n', [fam(:,1)*km, fam(:,4)*day, fam(:,5)].');
fprintf('\nTarget orbit (C = %.2f): Ax = %.0f km, x0 = %.10f, vy0 = %.10f\n', ...
        Ct, orbit(1)*km, orbit(2), orbit(3));
fprintf('Period = %.6f = %.3f days, |C - Ct| = %.1e\n', orbit(4), orbit(4)*day, abs(orbit(5) - Ct));
fprintf('Monodromy eigenvalues:\n'); fprintf('   %12.6g %+12.6gi\n', [real(lm), imag(lm)].');

% Save the target orbit for step 4c-2

if ~isfolder(fullfile(root, 'data')), mkdir(fullfile(root, 'data')); end

save(fullfile(root, 'data', 'lyapunovL1_C318.mat'), 'orbit', 'M', 'mu', 'Ct');

% ---------------------------------------------------------------------

fig = figure('Color', T.bg, 'Position', [100 100 1250 560]);
tl  = tiledlayout(fig, 1, 2, 'Padding', 'loose', 'TileSpacing', 'compact');

% (a) family over the zero-velocity curve of the target energy

ax1 = nexttile(tl); hold(ax1, 'on')
xv = linspace(0.76, 0.92, 641);  yv = linspace(-0.08, 0.08, 641);
[X, Y] = meshgrid(xv, yv);
TwoU = reshape(jacobiConstant([X(:), Y(:), zeros(numel(X), 4)], mu), size(X));
zvcPanel(ax1, T, xv, yv, X, Y, TwoU, Ct, mu, L)
nF  = size(fam, 1);
cols = [linspace(T.accent1(1), T.accent2(1), nF); ...
        linspace(T.accent1(2), T.accent2(2), nF); ...
        linspace(T.accent1(3), T.accent2(3), nF)].';   % small (blue) -> large (orange)
for k = 1:nF

    [~, S] = ode113(@(t,s) cr3bpEOM(t, s, mu), [0 fam(k,4)], [fam(k,2); 0; 0; 0; fam(k,3); 0], opts);
    plot(ax1, S(:,1), S(:,2), 'Color', cols(k,:), 'LineWidth', 0.8)

end

[~, S] = ode113(@(t,s) cr3bpEOM(t, s, mu), [0 orbit(4)], s0, opts);
hT = plot(ax1, S(:,1), S(:,2), 'Color', T.accent3, 'LineWidth', 2.2);
xlim(ax1, [0.76 0.92]); ylim(ax1, [-0.08 0.08])
legend(ax1, hT, sprintf('C = %.2f (step 3 energy)', Ct), 'TextColor', T.fg, ...
       'Color', T.bg, 'EdgeColor', T.grid, 'Location', 'southwest')
title(ax1, 'Lyapunov family around L1')

% (b) period against Jacobi constant

ax2 = nexttile(tl); hold(ax2, 'on')
plot(ax2, fam(:,5), fam(:,4)*day, 'o-', 'Color', T.accent2, 'MarkerFaceColor', T.accent2)
plot(ax2, orbit(5), orbit(4)*day, 'p', 'Color', T.accent3, 'MarkerFaceColor', T.accent3, 'MarkerSize', 12)
yline(ax2, 2*pi/nu*day, '--', 'linear limit 2\pi/\nu', 'Color', T.grid, 'LabelHorizontalAlignment', 'left')
xline(ax2, CL(1), ':', 'C_{L1}', 'Color', T.grid)
styleAxes(ax2, T); set(ax2, 'XDir', 'reverse')
xlabel(ax2, 'Jacobi constant C  (energy increases \rightarrow)'); ylabel(ax2, 'Period [days]')
title(ax2, 'Period along the family')

title(tl, 'Step 4c-1: Lyapunov family by continuation', 'Color', T.fg, 'FontName', T.font)
set(fig, 'InvertHardcopy', 'off')
print(fig, fullfile(root, 'figures', 'step4c_family.png'), '-dpng', '-r200')