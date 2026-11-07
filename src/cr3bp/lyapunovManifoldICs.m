function [Xmoon, Xearth, tOrb] = lyapunovManifoldICs(s0, Tp, mu, kind, N, epsilon)
%LYAPUNOVMANIFOLDICS  Initial states on the stable or unstable manifold of a periodic orbit.
%
%   [Xmoon, Xearth, tOrb] = lyapunovManifoldICs(s0, Tp, mu, kind, N, epsilon)
%
%   s0      : 6x1 initial state of the periodic orbit
%   Tp      : period of the orbit
%   kind    : 'unstable' or 'stable'
%   N       : number of points along the orbit
%   epsilon : displacement along the eigenvector (default 1e-6)
%   Xmoon   : N x 6 states displaced towards the Moon side (+x)
%   Xearth  : N x 6 states displaced towards the Earth side (-x)
%   tOrb    : N x 1 times along the orbit of each point
%
%   Unstable: integrate the states forwards. Stable: integrate backwards.
%
%   References:
%     Koon, Lo, Marsden & Ross, Dynamical Systems, the Three-Body Problem
%       and Space Mission Design, 2nd ed., 2022, Sec. 4.4, Step 4,
%       Eqs. (4.4.2)-(4.4.3), pp. 114-116.

    if nargin < 6, epsilon = 1e-6; end

    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    % Orbit + STM at N+1 equally spaced times (the last one is t = Tp)

    tOrb = linspace(0, Tp, N+1).';
    [~, W] = ode113(@(t,w) cr3bpEOMstm(t, w, mu), tOrb, [s0(:); reshape(eye(6), 36, 1)], opts);

    % Eigenvector of the monodromy matrix M = Phi(Tp)

    M = reshape(W(end, 7:42), 6, 6);
    [V, D] = eig(M);  lm = diag(D);

    if strcmpi(kind, 'unstable'), [~, i] = max(abs(lm)); else, [~, i] = min(abs(lm)); end

    Y0 = real(V(:, i));  Y0 = Y0/norm(Y0);

    if Y0(1) < 0, Y0 = -Y0; end              % convention: + side points to the Moon (+x)

    % Transport the eigenvector along the orbit and displace the states

    Xmoon = zeros(N, 6);  Xearth = zeros(N, 6);

    for k = 1:N

        s   = W(k, 1:6).';
        Phi = reshape(W(k, 7:42), 6, 6);
        Y   = Phi*Y0;  Y = Y/norm(Y);          % Eq. (4.4.3)
        Xmoon(k,:)  = (s + epsilon*Y).';       % Eq. (4.4.2)
        Xearth(k,:) = (s - epsilon*Y).';

    end

    tOrb = tOrb(1:N);
    
end