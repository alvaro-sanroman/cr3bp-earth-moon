function A = cr3bpJacobian(s, mu)

%CR3BPJACOBIAN  Jacobian of the CR3BP equations of motion, d(ds)/ds (6x6).
%
%   A = cr3bpJacobian(s, mu)
%
%   s  : state [x y z vx vy vz]
%   A  : [0 I; Uij 2*Omega], with Uij the second derivatives of U and
%        Omega = [0 1 0; -1 0 0; 0 0 0] (Coriolis terms of Eq. 3.50)
%
%   Second derivatives obtained by differentiating Wakker Eq. (3.49).
%   At a collinear point they reduce to Uxx = a, Uyy = -b (KLMR Eq. 2.6.2).
%
%   References:
%     K.F. Wakker, Fundamentals of Astrodynamics, TU Delft, 2015,
%       Sec. 3.3, Eqs. (3.49)-(3.50) (course reader).
%     W.S. Koon, M.W. Lo, J.E. Marsden, S.D. Ross, Dynamical Systems,
%       the Three-Body Problem and Space Mission Design, 2nd ed., 2022,
%       Secs. 2.6-2.7.

    x = s(1); y = s(2); z = s(3);
    dx1 = x + mu;  dx2 = x - 1 + mu;           % x-distance to each primary
    r1 = sqrt(dx1^2 + y^2 + z^2);
    r2 = sqrt(dx2^2 + y^2 + z^2);
    k1 = (1 - mu)/r1^3;  k2 = mu/r2^3;         % 1/r^3 terms
    c1 = 3*(1 - mu)/r1^5; c2 = 3*mu/r2^5;      % 3/r^5 terms

    Uxx = 1 - k1 - k2 + c1*dx1^2 + c2*dx2^2;
    Uyy = 1 - k1 - k2 + c1*y^2   + c2*y^2;
    Uzz =   - k1 - k2 + c1*z^2   + c2*z^2;
    Uxy = c1*dx1*y + c2*dx2*y;
    Uxz = c1*dx1*z + c2*dx2*z;
    Uyz = c1*y*z   + c2*y*z;

    Uij   = [Uxx Uxy Uxz; Uxy Uyy Uyz; Uxz Uyz Uzz];
    Omega = [0 1 0; -1 0 0; 0 0 0];
    A = [zeros(3), eye(3); Uij, 2*Omega];
    
end