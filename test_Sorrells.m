%% Propagator Matrix Method: P-SV system

% Benchmark with Sorrells solution for elastic halfspace

addpath('./src', './data');

%% Pressure source & Elastic properties

%%% Create pressure source %%%

% Mesh information [km]
mesh.Nx = 128;  mesh.dx = 0.04;
mesh.Ny = 128;  mesh.dy = 0.04;

% Source: Pressure wave (single Fourier mode)
% Unit: Pa, positive for downward normal traction
param.amp = 1;  param.Nw_x = 10;  param.Nw_y = 12;  param.fw = 0.1;
src = create_psfc(2, mesh, param);

% Wavelength of pressure wave
% Note: Sorrells solution at quasi-static limit only depends on wavelength
fprintf('Wavelength: %g km\n', src.wavelen);

% Only take one snapshot for benchmark
src.pp = squeeze(src.pp(:,:,1));  src.time = src.time(1);


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

%% Analytical Sorrells solution

% Depth query points [km]
mesh.Nz = 100;  zq = linspace(0, h0, mesh.Nz);  mesh.dz = zq(2) - zq(1);

% Sorrells solution
sol_sorl = calc_sorrells(src, zq, hs_prop);

%% Seismic modeling (Layered medium)

% Depth query points
src.depth_query = 1;  src.zq = zq;

% Only output displacement fields
src.include_stress = 0;

% Quasi-static modeling result
% (Now positive upward for vertical displ. to compare with Boussinesq solution)
[sol_pm, ~] = qs_model(src, elast_prop);  sol_pm.uz = -sol_pm.uz;

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
