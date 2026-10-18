function zvcPanel(ax, T, xv, yv, X, Y, TwoU, C, mu, L)
%ZVCPANEL  Draw the zero-velocity curve 2U = C in the plane z = 0.
%   Shades the forbidden region 2U < C and marks primaries and Li.
%   TwoU : 2U on the grid [X, Y] = meshgrid(xv, yv)
%   Reference: Wakker, Fundamentals of Astrodynamics, TU Delft, 2015,
%     Sec. 3.6, Eq. (3.59).

    forbid = TwoU < C;                                   % V^2 < 0 here
    fillC  = repmat(reshape(T.body, 1, 1, 3), size(forbid));
    image(ax, xv, yv, fillC, 'AlphaData', 0.9*forbid)
    contour(ax, X, Y, TwoU, [C C], 'Color', T.fg, 'LineWidth', 0.8)
    plot(ax, -mu,  0, 'o', 'MarkerSize', 8, 'MarkerFaceColor', T.accent1, 'MarkerEdgeColor', 'none')
    plot(ax, 1-mu, 0, 'o', 'MarkerSize', 4, 'MarkerFaceColor', T.fg,      'MarkerEdgeColor', 'none')
    plot(ax, L(:,1), L(:,2), '+', 'Color', T.accent3, 'MarkerSize', 7, 'LineWidth', 1.2)
    text(ax, L(:,1), L(:,2) + 0.13, {'L1','L2','L3','L4','L5'}, ...
         'Color', T.fg, 'FontName', T.font, 'FontSize', 8, ...
         'HorizontalAlignment', 'center')
    styleAxes(ax, T)
    set(ax, 'YDir', 'normal', 'Layer', 'top')
    axis(ax, 'equal'); xlim(ax, [-1.5 1.5]); ylim(ax, [-1.5 1.5])
    xlabel(ax, 'x'); ylabel(ax, 'y')
    
end