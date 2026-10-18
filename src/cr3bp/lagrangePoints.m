function [L, CL] = lagrangePoints(mu)

%LAGRANGEPOINTS  Positions and Jacobi constants of the five libration points.
%
%   [L, CL] = lagrangePoints(mu)
%
%   mu : mass parameter m2/(m1+m2)
%   L  : 5x2 matrix, row i = [x y] of point Li (all have z = 0)
%   CL : 5x1 vector, Jacobi constant at each point (body at rest)
%
%   L1, L2, L3: roots of dU/dx = 0 on the x-axis (y = z = 0), solved
%   numerically with fzero, one root per interval.
%   L4, L5: equilateral points, r1 = r2 = 1, x = 1/2 - mu, y = +-sqrt(3)/2.
%
%   References (course readers, Fundamentals of Astrodynamics, TU Delft):
%     K.F. Wakker, Fundamentals of Astrodynamics, TU Delft, 2015,
%       Sec. 3.7, Eqs. (3.62)-(3.68).

% dU/dx on the x-axis (Wakker Eq. 3.64, lecture 6j slide 9)

f = @(x) x - (1 - mu)*(mu + x)./abs(mu + x).^3 ...
               + mu*(1 - mu - x)./abs(1 - mu - x).^3;

d = 1e-6;                                % stay clear of the primaries
xL1 = fzero(f, [-mu + d,     1 - mu - d]);   % between Earth and Moon
xL2 = fzero(f, [ 1 - mu + d, 2         ]);   % beyond the Moon
xL3 = fzero(f, [-2,         -mu - d    ]);   % beyond the Earth

L = [xL1,      0;
     xL2,      0;
     xL3,      0;
     0.5 - mu,  sqrt(3)/2;                % L4 (leading the Moon)
     0.5 - mu, -sqrt(3)/2];               % L5 (trailing the Moon)

% Jacobi constant of a body at rest at each point (V = 0)

CL = jacobiConstant([L, zeros(5, 4)], mu);

end