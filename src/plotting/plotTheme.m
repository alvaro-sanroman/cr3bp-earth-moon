function T = plotTheme()
%PLOTTHEME  Colours used by every figure of the project (dark theme).
%   T = plotTheme();   then use T.bg, T.fg, T.accent1, ...

    T.bg      = [0.059 0.071 0.102];   % figure/axes background
    T.fg      = [0.862 0.882 0.910];   % text, axes lines
    T.grid    = [0.300 0.330 0.390];   % grid lines
    T.body    = [0.180 0.200 0.250];   % airfoil fill
    T.accent1 = [0.322 0.773 0.965];   % main line   (cyan)
    T.accent2 = [1.000 0.659 0.278];   % second line (amber)
    T.accent3 = [0.937 0.392 0.486];   % highlights  (coral)
    T.accent4 = [0.580 0.890 0.533];   % extra       (green)
    T.font    = 'Helvetica';
end