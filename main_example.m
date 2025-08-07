%% Propagator Matrix Method: P-SV system

% Example workflow for a layered medium

addpath('./my_func', './my_data');

%% Example pressure source

% Unit: Pa, positive for pushing downward

% Mesh information [km]
mesh.Nx = 256;  mesh.dx = 0.02;
mesh.Ny = 256;  mesh.dy = 0.02;

% Case 1: Gaussian load (param: amp, standard deviation, Gaussian center)
% Case 2: Single Fourier mode (param: amp, kx, ky)
% Case 3: Turbulent pressure field from Cloud Model 1
% Case 4: Delta load (white spectrum)
example_case = 2;

switch example_case

    case 1
        % Parameters (amp: Pa, sigma: km, x0/y0: km)
        param.amp = 1;  param.sigma_x = 0.01;  param.sigma_y = 0.01;
        param.x0 = mesh.Nx*mesh.dx/2;  param.y0 = mesh.Ny*mesh.dy/2;

    case 2
        % Parameters (amp: Pa, Nw: integer number of cycles)
        param.amp = 1;  param.Nw_x = 12;  param.Nw_y = 10;

    case {3, 4}
        param = [];
end

% Create input pressure source
src = create_psfc(example_case, mesh, param);

% Pressure source src is a struct with following fields:
%   xh: x-axis grid point (dim: Nx *  1)
%   yh: y-axis grid point (dim: 1  * Ny)
%   pp: surface pressure  (dim: Nx * Ny)

%% Example elastic structure

% Layer properties (rho, vp, vs, thickness)
% Units: g/cm^3, km/s, km/s, km

% Read from file
% elast_prop = readmatrix('./my_data/vel_model.csv', 'NumHeaderLines', 1);

% Manually create model
elast_prop = [1.6, 1.45, 0.27, 0.2; ...
              1.9, 1.9, 0.6, 0.2; ...
              2.0, 2.1, 0.8, 0.6; ...
              2.2, 2.4, 0.9, 0];

% Print layered model
Nlayer = size(elast_prop, 1) - 1;
row_names = cellstr(num2str((1:Nlayer)'))';
row_names{end+1} = 'Halfspace';
disp('Layered media:');
disp(array2table(elast_prop, ...
    "VariableNames", {'rho (g/cm^3)', 'Vp (km/s)', 'Vs (km/s)', 'Thickness (km)'}, ...
    "RowNames", row_names));

% Note: You may proportionally enlarge the grid size and the depth extent
% of your layered model, and still remain numerically stable.

% If the wavenumber k is large (e.g., focus on small scale loading), then
% the depth of the elastic halfspace should be small.

%% Seismic modeling (Layered medium)

% Depth query points [km]
mesh.Nz = 11;  mesh.dz = 0.1;  zq = (0:mesh.Nz-1)' .* mesh.dz;

% Solve displacement-stress vector toward the surface
ds_pm = solve_ds(src, elast_prop, zq);

% Evaluate coefficients at the surface
coeff_pm = solve_coeff(src, ds_pm);

% Numerical solution
sol_pm = calc_layer(ds_pm, coeff_pm, 'xyz', true);

%% 2D plot on a horizontal plane

% Depth index (ind_z = 1 for surface z = 0)
ind_z = 11;  fprintf('Depth: %g km\n', sol_pm.z(ind_z));

% Plot input pressure
% Positive for compression, opposite to stress convention!
load('rwb_cb.mat', 'mcolor');  ax_prop.cmap = flip(mcolor);
ax_prop.clabel = 'Pressure (Pa)';  ax_prop.title = 'Surface Pressure';
ax_prop.cmax = max(abs(src.pp), [], 'all');
plot_2d(src.xh, src.yh, src.pp', ax_prop);

% Plot displacement components
plot_result_2d(sol_pm, 'disp', ind_z);

% Plot stress components
plot_result_2d(sol_pm, 'stress', ind_z);