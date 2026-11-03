function [vy0, T, hist] = lyapunovCorrect(x0, vy0, mu, tol, maxIter)

%LYAPUNOVCORRECT  Planar Lyapunov orbit by differential correction.
%
%   [vy0, T, hist] = lyapunovCorrect(x0, vy0, mu, tol, maxIter)
%
%   Starts perpendicular to the x-axis, s0 = [x0 0 0 0 vy0 0], integrates
%   to the next x-axis crossing (half period) and corrects vy0 until vx = 0
%   there. A perpendicular crossing closes the orbit by the symmetry
%   y -> -y, t -> -t. x0 is kept fixed.
%
%   vy0  : corrected initial y-velocity
%   T    : period (twice the half-period crossing time)
%   hist : |vx| at the crossing, one value per iteration
%
%   References:
%     Koon, Lo, Marsden & Ross, Dynamical Systems, the Three-Body Problem
%       and Space Mission Design, 2nd ed., 2022, Sec. 4.4, Step 3 (pp. 111-114).

    if nargin < 4, tol = 1e-10;
    
    end
    
    if nargin < 5, maxIter = 20; 
    
    end
   
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @xAxisCrossing);
    hist = nan(maxIter, 1);

    for k = 1:maxIter
      
        w0 = [x0; 0; 0; 0; vy0; 0; reshape(eye(6), 36, 1)];
        [~, ~, te, we] = ode113(@(t,w) cr3bpEOMstm(t, w, mu), [0 10], w0, opts);

        t1  = te(end);                          % half period
        s1  = we(end, 1:6).';                   % state at the crossing
        Phi = reshape(we(end, 7:42), 6, 6);     % STM at the crossing
        vx1 = s1(4);  vy1 = s1(5);
        hist(k) = abs(vx1);

      
        if abs(vx1) < tol                       % perpendicular crossing: done
            
            T = 2*t1;  hist = hist(1:k);  return
       
        end

        ds1  = cr3bpEOM(0, s1, mu);  ax1 = ds1(4);          % dvx/dt at the crossing
        dvy0 = vx1 / (Phi(4,5) - Phi(2,5)*ax1/vy1);         % KLMR p. 113
        vy0  = vy0 - dvy0;
    
    end
   
    error('lyapunovCorrect:noConvergence', 'No convergence in %d iterations', maxIter)

end

% ---------------------------------------------------------------------

function [val, isterminal, direction] = xAxisCrossing(~, w)

%XAXISCROSSING  Event: y = 0 with y decreasing (half period of the orbit).
  
    val        = w(2);
    isterminal = 1;                             % stop the integration there
    direction  = -1;

end