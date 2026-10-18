function ds = cr3bpEOM(~, s, mu)

%CR3BPEOM  Equations of motion of the circular restricted three-body problem.
%
%   ds = cr3bpEOM(t, s, mu)
%
%   s  : state [x; y; z; vx; vy; vz] in the rotating, non-dimensional frame
%   mu : mass parameter m2/(m1+m2), mu <= 1/2  (Earth-Moon: ~0.01215)
%   ds : time derivative of the state [vx; vy; vz; ax; ay; az]
%
%   Frame: P1 (Earth) at x = -mu, P2 (Moon) at x = 1-mu, unit distance
%   P1P2, unit angular velocity. Time is not used (autonomous system).
%
%   References (course readers, Fundamentals of Astrodynamics, TU Delft):
%     K.F. Wakker, Fundamentals of Astrodynamics, TU Delft, 2015,
%       Sec. 3.3, Eqs. (3.48)-(3.50).
%     H.D. Curtis, Orbital Mechanics for Engineering Students, 3rd ed.,
%       Elsevier, 2014, Sec. 2.12.

x  = s(1);  y  = s(2);  z  = s(3);
vx = s(4);  vy = s(5);  vz = s(6);

% distances to P1 (Earth) and P2 (Moon)

r1 = sqrt((mu + x)^2     + y^2 + z^2);
r2 = sqrt((1 - mu - x)^2 + y^2 + z^2);

% partial derivatives of U = (x^2+y^2)/2 + (1-mu)/r1 + mu/r2, Eq. (3.49)
% (gravity of both primaries + centrifugal acceleration)

dUdx = x - (1 - mu)/r1^3*(mu + x) + mu/r2^3*(1 - mu - x);
dUdy = y - (1 - mu)/r1^3*y        - mu/r2^3*y;
dUdz =   - (1 - mu)/r1^3*z        - mu/r2^3*z;

% Eq. (3.50) solved for the accelerations (+-2v terms = Coriolis)

ax =  2*vy + dUdx;
ay = -2*vx + dUdy;
az =         dUdz;

ds = [vx; vy; vz; ax; ay; az];

end