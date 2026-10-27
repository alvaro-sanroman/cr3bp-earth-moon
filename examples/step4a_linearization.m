%% Step 4a - Linearization around L1 and L2: saddle x center x center
% References:
%   K.F. Wakker, Fundamentals of Astrodynamics, TU Delft, 2015,
%     Sec. 3.3, Eq. (3.49) (course reader).
%   Koon, Lo, Marsden & Ross, Dynamical Systems, the Three-Body Problem
%     and Space Mission Design, 2nd ed., 2022, Secs. 2.6-2.7,
%     Eqs. (2.6.2), (2.6.4) and p. 47.

clear; clc; close all
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')))

mu  = 1/(1 + 81.30056);
L   = lagrangePoints(mu);
day = 27.321661/(2*pi);                       % 1 time unit [days]
f   = @(s) cr3bpEOM(0, s, mu);

for i = 1:2
    
    sL = [L(i,1); 0; 0; 0; 0; 0];

    % Our Jacobian and its eigenvalues
    
    A  = cr3bpJacobian(sL, mu);
    ev = eig(A);
    lam = max(real(ev));                      % real pair +-lambda
    im  = sort(abs(imag(ev(abs(imag(ev)) > 1e-8))));   % +-nu_z, +-nu (each twice)
    nuZ = im(1);  nu = im(end);

    % KLMR closed form, Eq. (2.6.2) and p. 47
    
    mb  = mu/abs(L(i,1) - 1 + mu)^3 + (1 - mu)/abs(L(i,1) + mu)^3;
    a1  = (mb - 2 + sqrt(9*mb^2 - 8*mb))/2;
    a2  = (mb - 2 - sqrt(9*mb^2 - 8*mb))/2;
    lamK = sqrt(a1);  nuK = sqrt(-a2);  nuZK = sqrt(mb);

    % Finite-difference check of the Jacobian (central differences)
    
    h = 1e-6;  Afd = zeros(6);
    for j = 1:6
        e = zeros(6,1); e(j) = h;
        Afd(:,j) = (f(sL + e) - f(sL - e))/(2*h);
    end

    fprintf('L%d  (x = %.6f, mu_bar = %.4f)\n', i, L(i,1), mb);
    fprintf('  lambda = %.6f   KLMR %.6f   e-folding time 1/lambda = %.2f days\n', lam, lamK, day/lam);
    fprintf('  nu     = %.6f   KLMR %.6f   in-plane period 2pi/nu   = %.2f days\n', nu, nuK, 2*pi/nu*day);
    fprintf('  nu_z   = %.6f   sqrt(mu_bar) %.6f  out-of-plane period = %.2f days\n', nuZ, nuZK, 2*pi/nuZ*day);
    fprintf('  max |A - A_fd| = %.1e\n\n', max(abs(A - Afd), [], 'all'));

end