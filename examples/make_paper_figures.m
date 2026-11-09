function make_paper_figures()

%MAKE_PAPER_FIGURES  Re-run the example scripts with the light theme and
%   collect the figures in figures/paper/ for the AIAA paper.
%
%   The dark figures in figures/ (used in the README / GitHub) are left
%   exactly as they were, and so are the .mat files in data/.
%   Run it from anywhere:   >> make_paper_figures
%   Takes a few minutes (it runs every step again).

    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(root, 'src')))
    figDir   = fullfile(root, 'figures');
    paperDir = fullfile(figDir, 'paper');
    dataDir  = fullfile(root, 'data');

    if ~exist(paperDir, 'dir'), mkdir(paperDir); end

    % Scripts that produce static PNG figures (animations and video are skipped)

    scripts = { 'step1_eom', 'step2_zeroVelocity', 'step2_hillSurface', ...
                'step3_transit', 'step4a_linearization', 'step4a_mcgehee', ...
                'step4b_lyapunov', 'step4c_family', 'step4c_manifolds', ...
                'step4c_transitCheck' };

    % 1. Back up the dark figures and the data files

    bak = tempname;  mkdir(bak);  mkdir(fullfile(bak, 'fig'));  mkdir(fullfile(bak, 'data'));
    copyIfAny(fullfile(figDir, '*.png'), fullfile(bak, 'fig'));
    copyIfAny(fullfile(dataDir, '*.mat'), fullfile(bak, 'data'));

    % 2. Run every script with the light theme

    setappdata(groot, 'cr3bpTheme', 'light');
    restoreTheme = onCleanup(@() setappdata(groot, 'cr3bpTheme', 'dark'));
    made = {};

    for k = 1:numel(scripts)

        file = fullfile(root, 'examples', [scripts{k} '.m']);

        if ~isfile(file)

            fprintf('  (skipped, not found) %s\n', scripts{k});

            continue

        end

        fprintf('\n==== %s (light theme) ====\n', scripts{k});
        t0 = now;
        runOne(file);
        close all
        % move every PNG written by this script into figures/paper/

        d = dir(fullfile(figDir, '*.png'));

        for j = 1:numel(d)

            if d(j).datenum >= t0 - 2/86400

                movefile(fullfile(figDir, d(j).name), fullfile(paperDir, d(j).name), 'f');
                made{end + 1} = d(j).name; %#ok<AGROW>

            end
        end
    end

    % 3. Restore the dark figures and the data files

    copyIfAny(fullfile(bak, 'fig', '*.png'), figDir);
    copyIfAny(fullfile(bak, 'data', '*.mat'), dataDir);
    rmdir(bak, 's');
    clear restoreTheme                       % back to the dark theme

    fprintf('\nLight figures saved in %s:\n', paperDir);
    fprintf('  %s\n', made{:});

end

function runOne(file)

    % Own workspace: the 'clear' at the top of each script cannot touch
    % the variables of make_paper_figures.
    run(file);

end

function copyIfAny(pattern, dest)

    if ~isempty(dir(pattern))

        copyfile(pattern, dest, 'f');
        
    end
end