%% Function: Make comparison 2D pcolor plots

function plot_compare_2d(sol_bsnq, sol_pm, comp, ind_z)

    % Check spatial axis
    if any(sol_bsnq.x ~= sol_pm.x) || any(sol_bsnq.y ~= sol_pm.y)
        error('Spatial coordinates do not match!');
    end

    if strcmp(comp, 'disp')
        plot_compare_disp(sol_bsnq, sol_pm, ind_z);

    elseif strcmp(comp, 'stress')
        plot_compare_stress(sol_bsnq, sol_pm, ind_z);
    
    else
        error('Invalid component. Should be disp or stress!');
    end

end

%% Function: Compare displacement

function plot_compare_disp(sol_bsnq, sol_pm, iz)

    % Load colorbar (red-white-blue)
    load('rwb_cb.mat', 'mcolor');

    % Get the screen size
    screen = get(0, 'ScreenSize');

    % Color limit clipping factor
    cb_clip = 5;

    % Axes properties
    x = sol_bsnq.x;  y = sol_bsnq.y;
    ax_prop.cmap = mcolor;
    ax_prop.clabel = 'Displacement (μm)';
    
    % Make plots
    figure('Name', 'Displacement', 'Position', [0, 0, screen(3)*3/4, screen(4)]);
    tiledlayout(3, 3, 'Padding', 'loose', 'TileSpacing', 'loose');
    colormap(mcolor);
    
    %%% Displacement components %%%
    nexttile(1);  plot_var = sol_bsnq.uz(:,:,iz);
    ax_prop.title = 'Vertical';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(4);  plot_var = sol_pm.uz(:,:,iz);
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(2);  plot_var = sol_bsnq.ux(:,:,iz);
    ax_prop.title = 'X-dir';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(5);  plot_var = sol_pm.ux(:,:,iz);
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(3);  plot_var = sol_bsnq.uy(:,:,iz);
    ax_prop.title = 'Y-dir';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(6);  plot_var = sol_pm.uy(:,:,iz);
    plot_2d(x, y, plot_var, ax_prop);
    
    %%% Difference of displacement %%%
    nexttile(7);  plot_var = sol_bsnq.uz(:,:,iz) - sol_pm.uz(:,:,iz);
    ax_prop.title = 'Diff. (Vertical)';
    ax_prop.cmax = max(abs(plot_var), [], 'all');
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(8);  plot_var = sol_bsnq.ux(:,:,iz) - sol_pm.ux(:,:,iz);
    ax_prop.title = 'Diff. (X-dir)';
    ax_prop.cmax = max(abs(plot_var), [], 'all');
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(9);  plot_var = sol_bsnq.uy(:,:,iz) - sol_pm.uy(:,:,iz);
    ax_prop.title = 'Diff. (Y-dir)';
    ax_prop.cmax = max(abs(plot_var), [], 'all');
    plot_2d(x, y, plot_var, ax_prop);

end

%% Function: Compare stress

function plot_compare_stress(sol_bsnq, sol_pm, iz)

    % Load colorbar (red-white-blue)
    load('rwb_cb.mat', 'mcolor');  colormap(mcolor);

    % Get the screen size
    screen = get(0, 'ScreenSize');

    % Color limit clipping factor
    cb_clip = 25;

    % Axes properties
    x = sol_bsnq.x;  y = sol_bsnq.y;
    ax_prop.cmap = mcolor;
    ax_prop.clabel = 'Stress (Pa)';
    
    % Make plots
    figure('Name', 'Normal Stress', 'Position', [0, 0, screen(3)*3/4, screen(4)]);
    tiledlayout(3, 3, 'Padding', 'loose', 'TileSpacing', 'loose');
    colormap(mcolor);
    
    %%% Normal stress components %%%
    nexttile(1);  plot_var = sol_bsnq.szz(:,:,iz);
    ax_prop.title = 'Szz';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(4);  plot_var = sol_pm.szz(:,:,iz);
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(2);  plot_var = sol_bsnq.sxx(:,:,iz);
    ax_prop.title = 'Sxx';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(5);  plot_var = sol_pm.sxx(:,:,iz);
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(3);  plot_var = sol_bsnq.syy(:,:,iz);
    ax_prop.title = 'Syy';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(6);  plot_var = sol_pm.syy(:,:,iz);
    plot_2d(x, y, plot_var, ax_prop);
    
    %%% Difference of normal stress %%%
    nexttile(7);  plot_var = sol_bsnq.szz(:,:,iz) - sol_pm.szz(:,:,iz);
    ax_prop.title = 'Diff. (Szz)';
    ax_prop.cmax = max(abs(plot_var), [], 'all');
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(8);  plot_var = sol_bsnq.sxx(:,:,iz) - sol_pm.sxx(:,:,iz);
    ax_prop.title = 'Diff. (Sxx)';
    ax_prop.cmax = max(abs(plot_var), [], 'all');
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(9);  plot_var = sol_bsnq.syy(:,:,iz) - sol_pm.syy(:,:,iz);
    ax_prop.title = 'Diff. (Syy)';
    ax_prop.cmax = max(abs(plot_var), [], 'all');
    plot_2d(x, y, plot_var, ax_prop);
    
    % Make plots
    figure('Name', 'Shear Stress', 'Position', [0, 0, screen(3)*3/4, screen(4)]);
    tiledlayout(3, 3, 'Padding', 'loose', 'TileSpacing', 'loose');
    colormap(mcolor);
    
    %%% Shear stress components %%%
    nexttile(1);  plot_var = sol_bsnq.sxz(:,:,iz);
    ax_prop.title = 'Sxz';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(4);  plot_var = sol_pm.sxz(:,:,iz);
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(2);  plot_var = sol_bsnq.syz(:,:,iz);
    ax_prop.title = 'Syz';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(5);  plot_var = sol_pm.syz(:,:,iz);
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(3);  plot_var = sol_bsnq.sxy(:,:,iz);
    ax_prop.title = 'Sxy';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(6);  plot_var = sol_pm.sxy(:,:,iz);
    plot_2d(x, y, plot_var, ax_prop);
    
    %%% Difference of stress %%%
    nexttile(7);  plot_var = sol_bsnq.sxz(:,:,iz) - sol_pm.sxz(:,:,iz);
    ax_prop.title = 'Diff. (Sxz)';
    ax_prop.cmax = max(abs(plot_var), [], 'all');
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(8);  plot_var = sol_bsnq.syz(:,:,iz) - sol_pm.syz(:,:,iz);
    ax_prop.title = 'Diff. (Syz)';
    ax_prop.cmax = max(abs(plot_var), [], 'all');
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(9);  plot_var = sol_bsnq.sxy(:,:,iz) - sol_pm.sxy(:,:,iz);
    ax_prop.title = 'Diff. (Sxy)';
    ax_prop.cmax = max(abs(plot_var), [], 'all');
    plot_2d(x, y, plot_var, ax_prop);

end
