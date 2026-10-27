%% Step 4a (b) - McGehee representation of the equilibrium region near L1
% Phase space near L1 on a fixed energy surface, drawn as a spherical shell:
% inner sphere n1 (Earth side), outer sphere n2 (Moon side), Lyapunov orbit
% at the equator, stable/unstable tubes, transit and non-transit orbits.
% Built from the linear solution in eigenvector coordinates (xi, eta, zeta).
% References:
%   Koon, Lo, Marsden & Ross, Dynamical Systems, the Three-Body Problem
%     and Space Mission Design, 2nd ed., 2022, Sec. 2.7 (Eqs. 2.7.3-2.7.5,
%     Fig. 2.7.2), Sec. 2.8 (Eq. 2.8.3, transit time p. 53).
%   McGehee, R., Some homoclinic orbits for the restricted three-body
%     problem, PhD thesis, Univ. Wisconsin, 1969 (as presented in KLMR).

clear; clc; close all
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')))

mu = 1/(1 + 81.30056);
L  = lagrangePoints(mu);
T  = plotTheme();

% Saddle and center rates at L1 from our Jacobian (planar block)

A   = cr3bpJacobian([L(1,1); 0; 0; 0; 0; 0], mu);
ev  = eig(A([1 2 4 5], [1 2 4 5]));          % planar: x, y, vx, vy
lam = max(real(ev));
nu  = max(abs(imag(ev)));

% Region R: energy El = eps, |eta - xi| <= c   (KLMR Eqs. 2.7.4, 2.7.6)

eps = 1;  c = 2;
rL  = sqrt(2*eps/nu);                         % |zeta| on orbit and tubes, rho* = 2eps/nu
d0  = 0.05;                                   % start of the tubes, close to the orbit

fig = figure('Color', T.bg, 'Position', [100 100 1300 620]);
tl  = tiledlayout(fig, 1, 2, 'Padding', 'loose', 'TileSpacing', 'compact');
views = [-35 20; 0 0];

for k = 1:2
    
    ax = nexttile(tl); hold(ax, 'on')

    % Bounding spheres n1 (radius 1) and n2 (radius 2)
    
    [Xs, Ys, Zs] = sphere(40);
    
    for R = [1 2]
        
        surf(ax, R*Xs, R*Ys, R*Zs, 'FaceColor', T.fg, 'FaceAlpha', 0.05, ...
             'EdgeColor', T.grid, 'EdgeAlpha', 0.15)
    end

    % Lyapunov orbit: xi = eta = 0
    
    t  = linspace(0, 2*pi/nu, 200);
    P  = mcg(0*t, 0*t, rL*exp(-1i*nu*t), lam, nu, eps, c);
    hL = plot3(ax, P(1,:), P(2,:), P(3,:), 'Color', T.fg, 'LineWidth', 3);

    % Tubes: unstable (eta = 0, xi grows), stable (xi = 0, eta decays)
    
    [PH, TT] = meshgrid(linspace(0, 2*pi, 41), linspace(0, log(c/d0)/lam, 60));
    g = @(v) reshape(v, size(TT));
    
    for sgn = [1 -1]
       
        % unstable tube (eta = 0): surface + every 5th orbit as a line
        
        Pu = mcg(sgn*d0*exp(lam*TT(:)).', 0*TT(:).', rL*exp(1i*(PH(:) - nu*TT(:))).', lam, nu, eps, c);
        Xu = g(Pu(1,:));  Yu = g(Pu(2,:));  Zu = g(Pu(3,:));
        hu = plot3(ax, Xu(:,1:5:end), Yu(:,1:5:end), Zu(:,1:5:end), 'Color', T.accent3, 'LineWidth', 0.8);
        hU = hu(1);
        plot3(ax, Xu(:,1:5:end), Yu(:,1:5:end), Zu(:,1:5:end), 'Color', T.accent3, 'LineWidth', 0.8)

        % stable tube (xi = 0): same
        
        Ps = mcg(0*TT(:).', sgn*d0*exp(lam*TT(:)).', rL*exp(1i*(PH(:) + nu*TT(:))).', lam, nu, eps, c);
        Xt = g(Ps(1,:));  Yt = g(Ps(2,:));  Zt = g(Ps(3,:));
        hs = plot3(ax, Xt(:,1:5:end), Yt(:,1:5:end), Zt(:,1:5:end), 'Color', T.accent4, 'LineWidth', 0.8);
        hS = hs(1);
        plot3(ax, Xt(:,1:5:end), Yt(:,1:5:end), Zt(:,1:5:end), 'Color', T.accent4, 'LineWidth', 0.8)
   
    end

    % Transit (xi*eta > 0) and non-transit (xi*eta < 0) orbits entering through n2
    
    for chi = [0.1, -0.1]
        
        eta0 = (c + sqrt(c^2 + 4*chi))/2;  xi0 = chi/eta0;    % eta - xi = c at t = 0
        Tr   = log(eta0/abs(xi0))/lam;                        % time in R, KLMR p. 53
        rz   = sqrt(2/nu*(eps - lam*chi));                    % |zeta| from Eq. (2.7.4)
       
        for ph = [0 pi]
            
            t = linspace(0, Tr, 200);
            P = mcg(xi0*exp(lam*t), eta0*exp(-lam*t), rz*exp(1i*(ph - nu*t)), lam, nu, eps, c);
            
            if chi > 0
                
                hT = plot3(ax, P(1,:), P(2,:), P(3,:), 'Color', T.accent2, 'LineWidth', 2);
            
            else
                
                hN = plot3(ax, P(1,:), P(2,:), P(3,:), 'Color', T.accent1, 'LineWidth', 2);
            
            end
        
        end
   
    end

    axis(ax, 'equal', 'off'); view(ax, views(k,:))

end

lg = legend(ax, [hL hU hS hT hN], {'Lyapunov orbit', 'Unstable tube', 'Stable tube', ...
            'Transit orbit', 'Non-transit orbit'}, 'TextColor', T.fg, 'Color', T.bg, ...
            'EdgeColor', T.grid, 'Location', 'eastoutside');
title(tl, 'McGehee representation of the L1 region (inner sphere: Earth side, outer: Moon side)', ...
      'Color', T.fg, 'FontName', T.font)
set(fig, 'InvertHardcopy', 'off')
print(fig, fullfile(root, 'figures', 'step4a_mcgehee.png'), '-dpng', '-r200')

% ---------------------------------------------------------------------
function P = mcg(xi, eta, zeta, lam, nu, eps, c)

%MCG  Map eigenvector coordinates (xi, eta, zeta) to the McGehee shell.
%   For fixed s = eta - xi the energy surface is a 2-sphere (KLMR p. 49):
%   (lam/4)(eta+xi)^2 + (nu/2)|zeta|^2 = eps + (lam/4)(eta-xi)^2.
%   Normalising gives a unit vector; s sets the radius (n1 -> 1, n2 -> 2).
    
    s  = eta - xi;   h = eta + xi;
    Rs = sqrt(eps + lam/4*s.^2);
    u  = [sqrt(nu/2)*real(zeta); sqrt(nu/2)*imag(zeta); sqrt(lam)/2*h] ./ Rs;
    P  = (1 + (s + c)/(2*c)) .* u;

end