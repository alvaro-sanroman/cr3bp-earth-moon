function C = jacobiConstant(S, mu)

%JACOBICONSTANT  Jacobi constant of the circular restricted three-body problem.
%
%   C = jacobiConstant(S, mu)
%
%   S  : states, one per row: [x y z vx vy vz]  (N x 6)
%   mu : mass parameter m2/(m1+m2)
%   C  : Jacobi constant of each state (N x 1)
%
%   C = 2U - V^2 = x^2 + y^2 + 2(1-mu)/r1 + 2mu/r2 - V^2
%   V is the velocity relative to the rotating frame. C is conserved along
%   any trajectory, so it is used to check the numerical integration.
%
%   References (course readers, Fundamentals of Astrodynamics, TU Delft):
%     K.F. Wakker, Fundamentals of Astrodynamics, TU Delft, 2015,
%       Sec. 3.4, Eq. (3.53); Sec. 3.6, Eq. (3.59).

x = S(:,1);  y = S(:,2);  z = S(:,3);
V2 = S(:,4).^2 + S(:,5).^2 + S(:,6).^2;     % V^2 (rotating frame)

r1 = sqrt((mu + x).^2     + y.^2 + z.^2);
r2 = sqrt((1 - mu - x).^2 + y.^2 + z.^2);

C = x.^2 + y.^2 + 2*(1 - mu)./r1 + 2*mu./r2 - V2;
end