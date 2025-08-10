%% Function: Make 2D pcolor plot

function plot_2d(x, y, data, ax_prop)

    % Color axis range
    cmax = ax_prop.cmax;
    if cmax == 0;  cmax = 0.1;  end

    % Make plot
    % Input data has size Nx * Ny
    pcolor(x, y, data');  shading interp;
    colormap(ax_prop.cmap);  cb = colorbar;  clim([-cmax, cmax]);
    xlabel('X (km)');   ylabel('Y (km)');
    ylabel(cb, ax_prop.clabel);  title(ax_prop.title);
    axis equal;  axis tight;

end
