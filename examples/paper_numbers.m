%% PAPER NUMBERS - every value quoted in the AIAA paper, in one run
% Prints (and saves to paper_numbers.txt in the repo root) the numbers
% used in the abstract, the text and the tables of the paper, in the
% order of the paper:
%   0. Units                         (Table 1)
%   1. Jacobi constant drift         (Sec. II.B)
%   2. Libration points              (Table 2, Sec. III)
%   3. Transit experiment, C = 3.18  (Sec. IV, abstract)
%   4. Linear dynamics at L1         (Sec. V)
%   5. Small Lyapunov orbit, Ax=1e-3 (Sec. VI.B-C)
%   6. Lyapunov orbit with C = 3.18  (Table 3, abstract)
%   7. Manifold tubes                (Sec. VII.A)
%   8. Transit window from the tube  (Sec. VII.B, abstract)
%   9. Example trajectory, 28 deg    (Sec. VII.B)
% Takes a few minutes (sections 3 and 8 integrate many trajectories).

clear; clc
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')))

outTxt = fullfile(root, 'paper_numbers.txt');
if isfile(outTxt), delete(outTxt); end
diary(outTxt)
cleanupDiary = onCleanup(@() diary('off'));

%% 0. Units and constants
mu   = 1/(1 + 81.30056);
LU   = 384400;                          % length unit [km]
TU   = 27.321661/(2*pi);                % time unit [days]
VU   = LU/(TU*86400);                   % velocity unit [km/s]
tol  = 1e-12;                           % RelTol = AbsTol used here
opts = odeset('RelTol', tol, 'AbsTol', tol);
eom  = @(t, s) cr3bpEOM(t, s, mu);
xE = -mu;  xM = 1 - mu;

fprintf('================ 0. UNITS (Table 1) ================\n');
fprintf('mu              = %.10f\n', mu);
fprintf('length unit     = %.0f km\n', LU);
fprintf('time unit       = %.4f days\n', TU);
fprintf('velocity unit   = %.4f km/s\n', VU);
fprintf('ode113 tolerances (RelTol = AbsTol) = %.0e\n', tol);

%% 1. Jacobi constant drift (step 1 trajectory)
s0 = [0.8; 0; 0; 0; 0.3; 0];
[t1, S1] = ode113(eom, [0 30], s0, opts);
C1 = jac(S1, mu);
fprintf('\n================ 1. JACOBI DRIFT (Sec. II.B) ================\n');
fprintf('start (0.8, 0) with vy = 0.3, %.0f time units = %.0f days\n', t1(end), t1(end)*TU);
fprintf('C0              = %.6f\n', C1(1));
fprintf('max |C - C0|    = %.2e   (%d steps)\n', max(abs(C1 - C1(1))), numel(t1) - 1);

%% 2. Libration points
[L, ~] = lagrangePoints(mu);
if size(L, 1) < size(L, 2), L = L.'; end
fprintf('\n================ 2. LIBRATION POINTS (Table 2) ================\n');
fprintf('Point        x            y           C          gamma\n');
gam = [xM - L(1,1), L(2,1) - xM, xE - L(3,1)];
maxAcc = 0;
for i = 1:5
    sL = [L(i,1); L(i,2); 0; 0; 0; 0];
    CL = jac(sL.', mu);
    a = cr3bpEOM(0, sL, mu);
    maxAcc = max(maxAcc, max(abs(a(4:6))));
    if i <= 3
        fprintf('L%d   %12.7f  %11.7f   %.6f   %.6f\n', i, L(i,1), L(i,2), CL, gam(i));
    else
        fprintf('L%d   %12.7f  %11.7f   %.6f\n', i, L(i,1), L(i,2), CL);
    end
end
fprintf('max acceleration at the points = %.1e\n', maxAcc);
fprintf('L1 - Moon distance  = %.0f km\n', gam(1)*LU);
fprintf('L2 - Moon distance  = %.0f km\n', gam(2)*LU);
gN = fzero(@(g) (1 - mu)/(1 - g)^2 - mu/g^2, [0.01 0.5]);
fprintf('neutral point gamma_N = %.5f  (%.0f km from the Moon)\n', gN, gN*LU);
xL1 = L(1,1);

%% 3. Transit experiment at C = 3.18
C  = 3.18;
x0 = 0.80;
V0 = sqrt(twoU(x0, 0, mu) - C);
fprintf('\n================ 3. TRANSIT EXPERIMENT (Sec. IV) ================\n');
fprintf('C = %.2f,  start (%.2f, 0)\n', C, x0);
fprintf('V0              = %.4f  (%.0f m/s)\n', V0, V0*VU*1000);

tEnd  = 26/TU;
thC   = -175:5:180;                         % 72 directions, 5 deg
crossT = nan(size(thC));  strict = false(size(thC));  rMin = inf(size(thC));
for k = 1:numel(thC)
    sk = launch(x0, V0, thC(k));
    % (a) coarse criterion: crosses x = xL1 moving towards the Moon
    oe = odeset(opts, 'Events', @(t, s) evCross(t, s, xL1, +1));
    sol = ode113(eom, [0 tEnd], sk, oe);
    if ~isempty(sol.xe)
        crossT(k) = sol.xe(1)*TU;
        tt = linspace(0, tEnd, 20000);
        ss = deval(ode113(eom, [0 tEnd], sk, opts), tt);
        rMin(k) = min(hypot(ss(1,:) - xM, ss(2,:)))*LU;
    end
    % (b) strict criterion: reaches x = 1-mu before x = -mu (the Earth)
    strict(k) = isStrictTransit(sk, tEnd, mu, opts);
end
isCoarse = ~isnan(crossT);
fprintf('5 deg steps, 26 days:\n');
fprintf('  coarse (crosses xL1): %d of %d, theta = %s deg\n', nnz(isCoarse), numel(thC), mat2str(thC(isCoarse)));
fprintf('  crossing times        %.2f to %.2f days\n', min(crossT(isCoarse)), max(crossT(isCoarse)));
fprintf('  closest Moon pass     %.0f km from its centre\n', min(rMin(isCoarse)));
fprintf('  strict (reaches 1-mu before the Earth): %d of %d, theta = %s deg\n', ...
    nnz(strict), numel(thC), mat2str(thC(strict)));

dTh   = 1;                                  % fine scan step [deg]
thF   = -180:dTh:179;
okF   = false(size(thF));
for k = 1:numel(thF)
    okF(k) = isStrictTransit(launch(x0, V0, thF(k)), tEnd, mu, opts);
end
% contiguous window containing 28 deg
i28 = find(thF == 28);
a = i28;  while a > 1 && okF(a - 1), a = a - 1; end
b = i28;  while b < numel(thF) && okF(b + 1), b = b + 1; end
fprintf('fine scan, %g deg steps, strict criterion, 26 days:\n', dTh);
fprintf('  window: %g to %g deg  (width %g deg = %.1f%% of 360)\n', ...
    thF(a), thF(b), thF(b) - thF(a), (thF(b) - thF(a))/360*100);
others = thF(okF);  others = others(others < thF(a) | others > thF(b));
fprintf('  other transits outside the window: %s\n', mat2str(others));

%% 4. Linear dynamics at L1
c2  = (1 - mu)/abs(xL1 + mu)^3 + mu/abs(xL1 - 1 + mu)^3;
lam = sqrt((c2 - 2 + sqrt(9*c2^2 - 8*c2))/2);
nu  = sqrt((2 - c2 + sqrt(9*c2^2 - 8*c2))/2);
nuz = sqrt(c2);
tau = -(nu^2 + 2*c2 + 1)/(2*nu);
fprintf('\n================ 4. LINEAR DYNAMICS AT L1 (Sec. V) ================\n');
fprintf('c2 = %.4f   lambda = %.4f   nu = %.4f   nu_z = %.4f   tau = %.4f\n', c2, lam, nu, nuz, tau);
fprintf('exp(lambda)     = %.1f per time unit (%.2f days)\n', exp(lam), TU);
fprintf('linear period 2pi/nu = %.3f days\n', 2*pi/nu*TU);

%% 5. Small Lyapunov orbit, Ax = 1e-3
Ax1  = 1e-3;
x01  = xL1 - Ax1;
vyS  = -Ax1*nu*tau;
[vyC, T1, hist] = lyapunovCorrect(x01, vyS, mu);
T1 = fullPeriod(x01, vyC, T1, mu, opts);
[M1, ~] = monodromy(x01, vyC, T1, mu, opts);
e1 = eig(M1);
fprintf('\n================ 5. SMALL ORBIT Ax = 1e-3 (Sec. VI.B-C) ================\n');
fprintf('Ax = %.0f km\n', Ax1*LU);
fprintf('seed vy0 = %.6f   corrected = %.6f   change = %.1f%%\n', vyS, vyC, (vyC - vyS)/vyS*100);
disp('correction history (as returned by lyapunovCorrect):'); disp(hist)
fprintf('period  = %.3f days   C = %.6f\n', T1*TU, jac([x01 0 0 0 vyC 0], mu));
printEig(e1)
fprintf('exp(lambda*T)   = %.0f\n', exp(lam*T1));
fprintf('nu_z*T = %.3f rad -> cos = %.4f, sin = %.4f\n', nuz*T1, cos(nuz*T1), sin(nuz*T1));

%% 6. Lyapunov orbit with C = 3.18 (saved by step4c_family)
S   = load(fullfile(root, 'data', 'lyapunovL1_C318.mat'));
Ax  = S.orbit(1);  xo = S.orbit(2);  vyo = S.orbit(3);  Co = S.orbit(5);
Tp  = fullPeriod(xo, vyo, S.orbit(4), mu, opts);
[M, ~] = monodromy(xo, vyo, Tp, mu, opts);
eM = eig(M);
[~, Sorb] = ode113(eom, linspace(0, Tp, 2000), [xo 0 0 0 vyo 0], opts);
fprintf('\n================ 6. ORBIT C = 3.18 (Table 3, abstract) ================\n');
fprintf('C               = %.10f\n', Co);
fprintf('Ax              = %.6f  (%.0f km)\n', Ax, Ax*LU);
fprintf('x0, vy0         = %.7f, %.7f\n', xo, vyo);
fprintf('x extent        = %.0f km to +%.0f km about L1\n', (min(Sorb(:,1)) - xL1)*LU, (max(Sorb(:,1)) - xL1)*LU);
fprintf('y extent        = +/- %.0f km\n', max(abs(Sorb(:,2)))*LU);
fprintf('period          = %.6f  (%.3f days)\n', Tp, Tp*TU);
printEig(eM)

%% 7. Manifold tubes: time to reach the sections
Nm = 40;
[Wm, Ys0, Yu0] = orbitSTM(xo, vyo, Tp, Nm, M, mu, opts);
fprintf('\n================ 7. MANIFOLD TUBES (Sec. VII.A) ================\n');
fprintf('epsilon = 1e-6, %d points on the orbit\n', Nm);
names = {'unstable, Moon side  -> x = 1-mu', 'unstable, Earth side -> x = 0.80', ...
         'stable,   Moon side  -> x = 1-mu', 'stable,   Earth side -> x = 0.80'};
Ysel = {Yu0, Yu0, Ys0, Ys0};  sgn = [1 -1 1 -1];  tdir = [1 1 -1 -1];  xsec = [xM 0.80 xM 0.80];
for c = 1:4
    tt = nan(Nm, 1);
    for k = 1:Nm
        X = branchStart(Wm(k, :), Ysel{c}, sgn(c));
        oe = odeset(opts, 'Events', @(t, s) evCross(t, s, xsec(c), 0));
        [~, ~, te] = ode113(eom, [0 tdir(c)*10], X, oe);
        if ~isempty(te), tt(k) = abs(te(1))*TU; end
    end
    fprintf('  %s: %.1f to %.1f days (%d of %d reached)\n', names{c}, min(tt), max(tt), nnz(~isnan(tt)), Nm);
end

%% 8. Transit window predicted by the stable tube (cut at x = 0.80)
Nc = 400;
Wc = orbitSTM(xo, vyo, Tp, Nc, M, mu, opts);
cut = nan(Nc, 3);                                % [y, vx, vy] at x = 0.80
for k = 1:Nc
    X  = branchStart(Wc(k, :), Ys0, -1);         % stable, Earth side
    oe = odeset(opts, 'Events', @(t, s) evCross(t, s, 0.80, 0));
    [~, ~, te, ye] = ode113(eom, [0 -10], X, oe);
    if ~isempty(te), cut(k, :) = ye(1, [2 4 5]); end
end
thEdge = [];
for k = 1:Nc
    j = mod(k, Nc) + 1;
    if cut(k,1)*cut(j,1) < 0                     % y changes sign: crosses y = 0
        f  = cut(k,1)/(cut(k,1) - cut(j,1));
        vx = cut(k,2) + f*(cut(j,2) - cut(k,2));
        vy = cut(k,3) + f*(cut(j,3) - cut(k,3));
        thEdge(end + 1) = atan2d(vy, vx);        %#ok<SAGROW>
    end
end
thEdge = sort(thEdge);
fprintf('\n================ 8. TRANSIT WINDOW FROM THE TUBE (Sec. VII.B) ================\n');
fprintf('%d points on the orbit, stable tube cut at x = 0.80\n', Nc);
fprintf('edges of the window: %s deg\n', mat2str(round(thEdge, 2)));
if numel(thEdge) == 2
    inTube = thC >= thEdge(1) & thC <= thEdge(2);
    fprintf('  width %.2f deg\n', diff(thEdge));
    fprintf('  5-deg directions inside the tube: %s\n', mat2str(thC(inTube)));
    fprintf('  agree with strict criterion: %d of %d directions\n', nnz(inTube == strict), numel(thC));
end

%% 9. Example trajectory at 28 deg (the one in the video)
s28 = launch(x0, V0, 28);
sol = ode113(eom, [0 25/TU], s28, opts);
tt  = linspace(0, 25/TU, 50000);
ss  = deval(sol, tt);
fprintf('\n================ 9. TRAJECTORY AT 28 DEG (Sec. VII.B) ================\n');
fprintf('crosses x_L1    at %.3f days\n', tt(find(ss(1,:) > xL1, 1))*TU);
fprintf('crosses x = 1-mu at %.3f days\n', tt(find(ss(1,:) > xM, 1))*TU);
fprintf('closest Moon pass = %.0f km from its centre\n', min(hypot(ss(1,:) - xM, ss(2,:)))*LU);
fprintf('stays with x > xL1 after crossing: %d\n', all(ss(1, find(ss(1,:) > xL1, 1):end) > xL1));

fprintf('\nSaved to %s\n', outTxt);
clear cleanupDiary

%% ---------------- Local functions ----------------
function v = twoU(x, y, mu)
    v = x.^2 + y.^2 + 2*(1 - mu)./hypot(x + mu, y) + 2*mu./hypot(x - 1 + mu, y);
end

function C = jac(S, mu)
    % Jacobi constant for states in rows: C = 2U - V^2 (planar + z)
    r1 = sqrt((S(:,1) + mu).^2 + S(:,2).^2 + S(:,3).^2);
    r2 = sqrt((S(:,1) - 1 + mu).^2 + S(:,2).^2 + S(:,3).^2);
    C = S(:,1).^2 + S(:,2).^2 + 2*(1 - mu)./r1 + 2*mu./r2 - sum(S(:,4:6).^2, 2);
end

function s = launch(x0, V0, th)
    s = [x0; 0; 0; V0*cosd(th); V0*sind(th); 0];
end

function ok = isStrictTransit(s, tEnd, mu, opts)
    oe = odeset(opts, 'Events', @(t, q) evStrict(t, q, mu));
    [~, ~, ~, ~, ie] = ode113(@(t, q) cr3bpEOM(t, q, mu), [0 tEnd], s, oe);
    ok = ~isempty(ie) && ie(1) == 1;
end

function [v, term, dirn] = evStrict(~, s, mu)
    v    = [s(1) - (1 - mu); s(1) + mu];   % Moon line, Earth line
    term = [1; 1];
    dirn = [1; -1];
end

function [v, term, dirn] = evCross(~, s, xc, d)
    v = s(1) - xc;  term = 1;  dirn = d;
end

function T = fullPeriod(x0, vy0, T, mu, opts)
    [~, S] = ode113(@(t, s) cr3bpEOM(t, s, mu), [0 T], [x0 0 0 0 vy0 0], opts);
    if norm(S(end, :) - [x0 0 0 0 vy0 0]) > 1e-6, T = 2*T; end
end

function [M, W] = monodromy(x0, vy0, T, mu, opts)
    w0 = [x0; 0; 0; 0; vy0; 0; reshape(eye(6), [], 1)];
    [~, W] = ode113(@(t, w) cr3bpEOMstm(t, w, mu), [0 T], w0, opts);
    M = reshape(W(end, 7:42), 6, 6);
end

function [W, Ys0, Yu0] = orbitSTM(x0, vy0, T, N, M, mu, opts)
    % State + STM at N points along the orbit, and the eigenvectors of M
    % normalised with positive x component ('+' = Moon side)
    w0 = [x0; 0; 0; 0; vy0; 0; reshape(eye(6), [], 1)];
    [~, W] = ode113(@(t, w) cr3bpEOMstm(t, w, mu), linspace(0, T, N + 1), w0, opts);
    W = W(1:N, :);
    [V, D] = eig(M);  l = diag(D);
    [~, iu] = max(abs(l));  [~, is] = min(abs(l));
    Yu0 = real(V(:, iu));  Yu0 = Yu0/norm(Yu0);  if Yu0(1) < 0, Yu0 = -Yu0; end
    Ys0 = real(V(:, is));  Ys0 = Ys0/norm(Ys0);  if Ys0(1) < 0, Ys0 = -Ys0; end
end

function X = branchStart(wRow, Y0, sgn)
    s   = wRow(1:6).';
    Phi = reshape(wRow(7:42), 6, 6);
    Y   = Phi*Y0;  Y = Y/norm(Y);
    X   = s + sgn*1e-6*Y;
end

function printEig(e)
    [~, i] = sort(abs(e), 'descend');  e = e(i);
    fprintf('eigenvalues of M:\n');
    fprintf('   %12.5f %+.5fi   |.| = %.6f\n', [real(e) imag(e) abs(e)].');
    r = e(abs(imag(e)) < 1e-9);
    fprintf('lambda1 = %.1f   lambda2 = %.3e   lambda1*lambda2 = %.6f\n', max(r), min(r), max(r)*min(r));
end