%% Propagator Matrix Method: P-SV system

% Benchmark with Boussinesq solution for elastic halfspace

addpath('./src', './data');

%% Pressure source & Elastic properties

%%% Create pressure source %%%

% Mesh information [km]
mesh.Nx = 256;  mesh.dx = 0.02;
mesh.Ny = 256;  mesh.dy = 0.02;

% Source: Point load from white spectrum
% Unit: Pa, positive for downward normal traction
src = create_psfc(4, mesh, []);

% Source: Gaussian load
% param.amp = 1;  param.sigma_x = 0.05;  param.sigma_y = 0.05;
% param.x0 = mesh.Nx*mesh.dx/2;  param.y0 = mesh.Ny*mesh.dy/2;
% src = create_psfc(1, mesh, param);

% Total force of point / Gaussian load
fprintf('Total force: %g N\n', sum(src.pp,'all')*mesh.dx*mesh.dy*1e6);


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

%% Analytical Boussinesq solution

% Depth query points [km]
mesh.Nz = 6;  mesh.dz = 0.025;  zq = (0:mesh.Nz-1)' .* mesh.dz;

% Remove query points below h0
zq = zq(zq <= h0);  mesh.Nz = length(zq);

% Boussinesq solution
sol_bsnq = calc_boussinesq(mesh, hs_prop, 'xyz');

%% Seismic modeling (Layered medium)

% Depth query points
src.depth_query = 1;  src.zq = zq;

% Output stress fields
src.include_stress = 1;

% Quasi-static modeling result
% (Now positive upward for vertical displ. to compare with Boussinesq solution)
[sol_pm, ~] = qs_model(src, elast_prop);  sol_pm.uz = -sol_pm.uz;

% Add an offset for comparison
% (the mean of Boussinesq solution in this domain)
varnames = {'ux', 'uy', 'uz', 'sxx', 'syy', 'szz', 'sxy', 'sxz', 'syz'};
for ivar = 1:numel(varnames)
   tmp1 = sol_pm.(varnames{ivar});  tmp2 = sol_bsnq.(varnames{ivar});
   sol_pm.(varnames{ivar}) = tmp1 + mean(tmp2-tmp1, [1,2], 'omitnan');
end
clear tmp1 tmp2;

%% Comparison on a horizontal plane: 2D plots

% Depth index (ind_z = 1 for surface z = 0)
ind_z = 3;  fprintf('Depth: %g km\n', sol_pm.zq(ind_z));

% Compare displacement components
plot_compare_2d(sol_bsnq, sol_pm, 'disp', ind_z);

% Compare stress components
plot_compare_2d(sol_bsnq, sol_pm, 'stress', ind_z);
