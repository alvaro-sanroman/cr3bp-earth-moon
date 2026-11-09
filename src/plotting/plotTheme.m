function T = plotTheme(mode)
%PLOTTHEME  Colours used by every figure of the project.
%   T = plotTheme()          current theme (dark unless switched to light)
%   T = plotTheme('dark')    dark theme  (GitHub / LinkedIn figures)
%   T = plotTheme('light')   light theme (paper figures)
%
%   The current theme is stored in the MATLAB session, so every script
%   that calls plotTheme() follows it without being edited:
%       setappdata(groot, 'cr3bpTheme', 'light')   % switch to light
%       setappdata(groot, 'cr3bpTheme', 'dark')    % back to dark
%   make_paper_figures.m does this automatically.

    if nargin < 1 || isempty(mode)

        mode = getappdata(groot, 'cr3bpTheme');
        
        if isempty(mode), mode = 'dark'; end

    end

    switch lower(mode)

        case 'light'

            T.bg      = [1.000 1.000 1.000];   % figure/axes background
            T.fg      = [0.100 0.100 0.120];   % text, axes lines
            T.grid    = [0.750 0.760 0.790];   % grid lines
            T.body    = [0.860 0.870 0.890];   % forbidden region, fills
            T.accent1 = [0.000 0.400 0.700];   % main line   (blue)
            T.accent2 = [0.850 0.400 0.000];   % second line (orange)
            T.accent3 = [0.780 0.150 0.300];   % highlights  (crimson)
            T.accent4 = [0.150 0.550 0.200];   % extra       (green)

        otherwise

            T.bg      = [0.059 0.071 0.102];
            T.fg      = [0.862 0.882 0.910];
            T.grid    = [0.300 0.330 0.390];
            T.body    = [0.180 0.200 0.250];
            T.accent1 = [0.322 0.773 0.965];   % cyan
            T.accent2 = [1.000 0.659 0.278];   % amber
            T.accent3 = [0.937 0.392 0.486];   % coral
            T.accent4 = [0.580 0.890 0.533];   % green
            mode = 'dark';

    end

    T.font = 'Helvetica';
    T.name = lower(mode);

end