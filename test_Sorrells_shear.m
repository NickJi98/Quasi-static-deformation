%% Propagator Matrix Method: P-SV + SH system

% Benchmark with Sorrells-type solution for shear traction wave

addpath('./src', './data');

%% Traction source & Elastic properties

%%% Create traction source %%%

% Mesh information [km]
mesh.Nx = 128;  mesh.dx = 0.04;
mesh.Ny = 128;  mesh.dy = 0.04;

% Source: Traction wave (single Fourier mode)
% Unit: Pa, positive along traction azimuth taz (here +x direction)
param.amp = 1;  param.Nw_x = 10;  param.Nw_y = 12;  param.fw = 0.1;
param.taz = 0;
src = create_tsfc(2, mesh, param);

% Wavelength of traction wave
% Note: Sorrells-type solution at static limit only depends on wavelength
fprintf('Wavelength: %g km\n', src.wavelen);
fprintf('Wave azimuth: %g deg, Traction azimuth: %g deg\n', ...
    rad2deg(src.waveaz), rad2deg(src.taz));

% Only take one snapshot for benchmark
src.tx = squeeze(src.tx(:,:,1));  src.ty = squeeze(src.ty(:,:,1));
src.time = src.time(1);


%%% Halfspace properties (rho, vp, vs) %%%
% Units: g/cm^3, km/s, km/s, km
hs_prop = [1.6, 1.45, 0.27];
% hs_prop = [2.7, 6.8, 4.0];

% Top layer thickness for propagator matrix method
% Unit: km, value between 0 to ~2 km
h0 = 0.6;
% Note: For halfspace benchmark, smaller h0 gives better accuracy
% Note: h0 should be greater than the depth to evaluate outputs

% Create artificial layered model for propagator matrix method
elast_prop = [hs_prop h0; hs_prop 0];
Nlayer = size(elast_prop, 1) - 1;

% Print layered model
row_names = cellstr(num2str((1:Nlayer)'))';
row_names{end+1} = 'Halfspace';
disp('Layered media:');
disp(array2table(elast_prop, ...
    "VariableNames", {'rho (g/cm^3)', 'Vp (km/s)', 'Vs (km/s)', 'Thickness (km)'}, ...
    "RowNames", row_names));

%% Analytical Sorrells-type solution

% Depth query points [km]
mesh.Nz = 100;  zq = linspace(0, h0, mesh.Nz);  mesh.dz = zq(2) - zq(1);

% Sorrells-type solution for shear traction
sol_sorl = calc_sorrells_shear(src, zq, hs_prop);

%% Seismic modeling (Layered medium)

% Depth query points
src.depth_query = 1;  src.zq = zq;

% Only output displacement fields
src.include_stress = 0;

% Quasi-static modeling result: P-SV & SH systems
% (Now positive upward for vertical displ. to compare with analytical solution)
[sol_psv, ~] = qs_model(src, elast_prop);  sol_psv.uz = -sol_psv.uz;
[sol_sh, ~] = qs_model_sh(src, elast_prop);

% Superpose P-SV and SH responses
varnames = {'ux', 'uy', 'uz'};
sol_pm = sol_psv;
for ivar = 1:numel(varnames)
   sol_pm.(varnames{ivar}) = sol_psv.(varnames{ivar}) + sol_sh.(varnames{ivar});
end

%% Plot depth profile

% Depth profile of amplitude
ux_max = squeeze(max(abs(sol_pm.ux), [], [1,2]));
uy_max = squeeze(max(abs(sol_pm.uy), [], [1,2]));
uz_max = squeeze(max(abs(sol_pm.uz), [], [1,2]));

% Plot depth profile
screen = get(0, 'ScreenSize');
figure('Name', 'Depth profile', 'Position', [0, 0, screen(3)/2, screen(4)/2.5]);
subplot(1,2,1);  plot(uz_max, zq, 'k-');  hold on;
plot(sol_sorl.uz, zq, 'r--');
xlabel('Displacement (μm)');  ylabel('Depth (km)');
set(gca, 'YDir', 'reverse');  grid off;
yline(src.wavelen, 'm-', 'LineWidth', 2);  title('Vertical');
legend('Numerical', 'Analytical', 'Location', 'best');

subplot(1,2,2);  plot(ux_max, src.zq, 'k-', uy_max, src.zq, 'b-');  hold on;
plot(sol_sorl.ux, zq, 'r--', sol_sorl.uy, zq, 'r-.');
xlabel('Displacement (μm)');  ylabel('Depth (km)');
set(gca, 'YDir', 'reverse');  grid off;
yline(src.wavelen, 'm-', 'LineWidth', 2);  title('Horizontal');
legend('Numerical (X)', 'Numerical (Y)', 'Analytical (X)', 'Analytical (Y)', ...
    'Location', 'best');

%% Quantitative comparison

% Maximum difference relative to peak amplitude (over depth profile)
fprintf('\nMax difference relative to peak amplitude (depth profile):\n');
fprintf('  ux : %.3e\n', max(abs(ux_max - sol_sorl.ux')) / max(sol_sorl.ux));
fprintf('  uy : %.3e\n', max(abs(uy_max - sol_sorl.uy')) / max(sol_sorl.uy));
fprintf('  uz : %.3e\n', max(abs(uz_max - sol_sorl.uz')) / max(sol_sorl.uz));

%% Save benchmark figures (./latex/Figures)

fig_dir = './latex/Figures';
if ~exist(fig_dir, 'dir');  mkdir(fig_dir);  end

fig = findobj('Type', 'figure', 'Name', 'Depth profile');
exportgraphics(fig(1), fullfile(fig_dir, 'Sorrells_shear.png'), 'Resolution', 150);
