%% Propagator Matrix Method: P-SV + SH system

% Benchmark with Cerruti solution for elastic halfspace

addpath('./src', './data');

%% Traction source & Elastic properties

%%% Create traction source %%%

% Mesh information [km]
mesh.Nx = 256;  mesh.dx = 0.02;
mesh.Ny = 256;  mesh.dy = 0.02;

% Source: Point load from white spectrum
% Unit: Pa, positive along traction azimuth taz (here +x direction)
param.taz = 0;
src = create_tsfc(4, mesh, param);

% Total force of point load
fprintf('Total force: %g N\n', sum(src.tx,'all')*mesh.dx*mesh.dy*1e6);


%%% Halfspace properties (rho, vp, vs) %%%
% Units: g/cm^3, km/s, km/s, km
hs_prop = [1.6, 1.45, 0.27];
% hs_prop = [2.7, 6.8, 4.0];

% Top layer thickness for propagator matrix method
% Unit: km, value between 0 to ~2 km
h0 = 0.1;
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

%% Analytical Cerruti solution

% Depth query points [km]
mesh.Nz = 6;  mesh.dz = 0.025;  zq = (0:mesh.Nz-1)' .* mesh.dz;

% Remove query points below h0
zq = zq(zq <= h0);  mesh.Nz = length(zq);

% Cerruti solution
sol_crt = calc_cerruti(mesh, hs_prop, 'xyz');

%% Seismic modeling (Layered medium)

% Depth query points
src.depth_query = 1;  src.zq = zq;

% Output stress fields
src.include_stress = 1;

% Quasi-static modeling result: P-SV & SH systems
% (Now positive upward for vertical displ. to compare with Cerruti solution)
[sol_psv, ~] = qs_model(src, elast_prop);  sol_psv.uz = -sol_psv.uz;
[sol_sh, ~] = qs_model_sh(src, elast_prop);

% Superpose P-SV and SH responses
varnames = {'ux', 'uy', 'uz', 'sxx', 'syy', 'szz', 'sxy', 'sxz', 'syz'};
sol_pm = sol_psv;
for ivar = 1:numel(varnames)
   sol_pm.(varnames{ivar}) = sol_psv.(varnames{ivar}) + sol_sh.(varnames{ivar});
end

% Add an offset for comparison
% (the mean of Cerruti solution in this domain)
for ivar = 1:numel(varnames)
   tmp1 = sol_pm.(varnames{ivar});  tmp2 = sol_crt.(varnames{ivar});
   sol_pm.(varnames{ivar}) = tmp1 + mean(tmp2-tmp1, [1,2], 'omitnan');
end
clear tmp1 tmp2;

%% Comparison on a horizontal plane: 2D plots

% Depth index (ind_z = 1 for surface z = 0)
ind_z = 3;  fprintf('Depth: %g km\n', sol_pm.zq(ind_z));

% Compare displacement components
plot_compare_2d(sol_crt, sol_pm, 'disp', ind_z);

% Compare stress components
plot_compare_2d(sol_crt, sol_pm, 'stress', ind_z);

%% Quantitative comparison

% Maximum difference relative to peak amplitude (at depth ind_z)
fprintf('\nMax difference relative to peak amplitude (depth %g km):\n', sol_pm.zq(ind_z));
for ivar = 1:numel(varnames)
    tmp1 = sol_pm.(varnames{ivar})(:,:,ind_z);
    tmp2 = sol_crt.(varnames{ivar})(:,:,ind_z);
    fprintf('  %3s : %.3e\n', varnames{ivar}, ...
        max(abs(tmp1-tmp2), [], 'all') / max(abs(tmp2), [], 'all'));
end
clear tmp1 tmp2;

%% Save benchmark figures (./latex/Figures)

fig_dir = './latex/Figures';
if ~exist(fig_dir, 'dir');  mkdir(fig_dir);  end

fig = findobj('Type', 'figure', 'Name', 'Displacement');
exportgraphics(fig(1), fullfile(fig_dir, 'Cerruti_disp.png'), 'Resolution', 150);
fig = findobj('Type', 'figure', 'Name', 'Normal Stress');
exportgraphics(fig(1), fullfile(fig_dir, 'Cerruti_normal_stress.png'), 'Resolution', 150);
fig = findobj('Type', 'figure', 'Name', 'Shear Stress');
exportgraphics(fig(1), fullfile(fig_dir, 'Cerruti_shear_stress.png'), 'Resolution', 150);
