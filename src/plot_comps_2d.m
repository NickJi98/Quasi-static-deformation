%% Function: 2D pcolor plots for numerical result

function plot_comps_2d(sol_pm, comp, inds)

    if strcmp(comp, 'disp')
        plot_result_disp(sol_pm, inds);

    elseif strcmp(comp, 'stress')
        plot_result_stress(sol_pm, inds);
    
    else
        error('Invalid component. Should be disp or stress!');
    end

end

%% Function: Compare displacement

function plot_result_disp(sol_pm, inds)

    % Load colorbar (red-white-blue)
    load('rwb_cb.mat', 'mcolor');

    % Get the screen size
    screen = get(0, 'ScreenSize');

    % Color limit clipping factor
    cb_clip = 1;

    % Depth and time indices
    iz = inds(1);  it = inds(2);

    % Axes properties
    x = sol_pm.x;  y = sol_pm.y;
    ax_prop.cmap = mcolor;
    ax_prop.clabel = 'Displacement (μm)';
    
    % Make plots
    figure('Name', 'Displacement', 'Position', [0, 0, screen(3)*3/4, screen(4)/3]);
    tiledlayout(1, 3, 'Padding', 'loose', 'TileSpacing', 'loose');
    colormap(mcolor);
    
    %%% Displacement components %%%
    nexttile(1);  plot_var = sol_pm.uz(:,:,iz,it);
    ax_prop.title = 'Vertical';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(2);  plot_var = sol_pm.ux(:,:,iz,it);
    ax_prop.title = 'X-dir';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(3);  plot_var = sol_pm.uy(:,:,iz,it);
    ax_prop.title = 'Y-dir';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);

end

%% Function: Compare stress

function plot_result_stress(sol_pm, inds)

    % Load colorbar (red-white-blue)
    load('rwb_cb.mat', 'mcolor');  colormap(mcolor);

    % Get the screen size
    screen = get(0, 'ScreenSize');

    % Color limit clipping factor
    cb_clip = 1;

    % Depth and time indices
    iz = inds(1);  it = inds(2);

    % Axes properties
    x = sol_pm.x;  y = sol_pm.y;
    ax_prop.cmap = mcolor;
    ax_prop.clabel = 'Stress (Pa)';
    
    % Make plots
    figure('Name', 'Normal Stress', 'Position', [0, screen(4)/3, screen(3)*3/4, screen(4)*2/3]);
    tiledlayout(2, 3, 'Padding', 'loose', 'TileSpacing', 'loose');
    colormap(mcolor);
    
    %%% Normal stress components %%%
    nexttile(1);  plot_var = sol_pm.szz(:,:,iz,it);
    ax_prop.title = 'Szz';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(2);  plot_var = sol_pm.sxx(:,:,iz,it);
    ax_prop.title = 'Sxx';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(3);  plot_var = sol_pm.syy(:,:,iz,it);
    ax_prop.title = 'Syy';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    %%% Shear stress components %%%
    nexttile(4);  plot_var = sol_pm.sxz(:,:,iz,it);
    ax_prop.title = 'Sxz';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(5);  plot_var = sol_pm.syz(:,:,iz,it);
    ax_prop.title = 'Syz';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);
    
    nexttile(6);  plot_var = sol_pm.sxy(:,:,iz,it);
    ax_prop.title = 'Sxy';
    ax_prop.cmax = max(abs(plot_var), [], 'all') / cb_clip;
    plot_2d(x, y, plot_var, ax_prop);

end
