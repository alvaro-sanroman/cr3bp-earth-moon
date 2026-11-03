function dw = cr3bpEOMstm(t, w, mu)

%CR3BPEOMSTM  CR3BP equations of motion plus variational equations (STM).
%
%   dw = cr3bpEOMstm(t, w, mu)
%
%   w : 42x1 = [s; Phi(:)], state (6) and state transition matrix (6x6,
%       stored column by column)
%   dPhi/dt = A(s) * Phi,  Phi(0) = I
%
%   References:
%     Koon, Lo, Marsden & Ross, Dynamical Systems, the Three-Body Problem
%       and Space Mission Design, 2nd ed., 2022, Sec. 4.2, Eqs. (4.2.3)-(4.2.5).

    s   = w(1:6);
    Phi = reshape(w(7:42), 6, 6);

    ds   = cr3bpEOM(t, s, mu);                % trajectory
    dPhi = cr3bpJacobian(s, mu) * Phi;        % variational equations, Eq. (4.2.5)

    dw = [ds; dPhi(:)];

end