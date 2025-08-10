%% Propagator Matrix Method: P-SV system

% Example workflow for a layered medium

% Add modeling functions
addpath('./src', './data');

% Parallel environment
% parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));

%% Example pressure source

% Unit: Pa, positive for downward normal traction

% Mesh information [km]
mesh.Nx = 128;  mesh.dx = 0.04;
mesh.Ny = 128;  mesh.dy = 0.04;

% Case 1: Gaussian load (param: amp, standard deviation, Gaussian center)
% Case 2: Single Fourier mode (param: amp, kx, ky, f)
% Case 3: Turbulent pressure field from Cloud Model 1
% Case 4: Delta load (white spectrum)
example_case = 2;

switch example_case

    case 1
        % Parameters (amp: Pa, sigma: km, x0/y0: km)
        param.amp = 1;  param.sigma_x = 0.01;  param.sigma_y = 0.01;
        param.x0 = mesh.Nx*mesh.dx/2;  param.y0 = mesh.Ny*mesh.dy/2;

    case 2
        % Parameters (amp: Pa, Nw: integer number of cycles, fw: Hz)
        param.amp = 1;
        param.Nw_x = 10;  param.Nw_y = 12;
        param.fw = 0.1;

    case {3, 4}
        param = [];
end

% Create input pressure source
src = create_psfc(example_case, mesh, param);

% Pressure source src is a struct with following fields:
%   xh: x-axis grid point (dim: Nx *  1)
%   yh: y-axis grid point (dim: 1  * Ny)
%   time: time samples    (dim: Nt, optional)
%   pp: surface pressure  (dim: Nx * Ny  or  Nx * Ny * Nt)

% Only Case 2 & 3 are time dependent

%% Example elastic structure

% Layer properties (rho, vp, vs, thickness)
% Units: g/cm^3, km/s, km/s, km

% Read from file
% elast_prop = readmatrix('./data/vel_model.csv', 'NumHeaderLines', 1);

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

% Optional: Depth query points [km]
src.depth_query = 1;
Nz = 5;  dz = 0.1;  src.zq = (0:Nz-1)' .* dz;

% Optional: Output stress fields
src.include_stress = 1;

% Quasi-static modeling result
[sol_pm, ~] = qs_model(src, elast_prop);

%% 2D plot on a horizontal plane

% Depth index (ind_z = 1 for surface z = 0)
ind_z = 2;
if length(sol_pm.zq) < ind_z
    ind_z = length(sol_pm.zq);
end
fprintf('Depth: %g km\n', sol_pm.zq(ind_z));

% Time index
ind_t = 3;
if isfield(sol_pm, 'time')
    fprintf('Time: %g s\n', sol_pm.time(ind_t));
else
    ind_t = 1;
    fprintf('Static modeling, set ind_t = 1\n');
end

% Plot input pressure
% Positive for compression, opposite to stress convention!
load('rwb_cb.mat', 'mcolor');  ax_prop.cmap = mcolor;
ax_prop.clabel = 'Pressure (Pa)';  ax_prop.title = 'Surface Pressure';
ax_prop.cmax = max(abs(src.pp), [], 'all');
plot_2d(src.xh, src.yh, src.pp(:,:,1), ax_prop);

% Plot displacement components
plot_comps_2d(sol_pm, 'disp', [ind_z, ind_t]);

% Plot stress components
plot_comps_2d(sol_pm, 'stress', [ind_z, ind_t]);
