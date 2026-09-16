%% Function: Solve displacement-stress vector (P-SV system)

% Return output vectors at different depths specified by zq
% Output at surface z = 0 is always recorded

function output = solve_ds(src, elast_prop)

    %%% Insert depth query points into layered model %%%
    model_prop = insert_query_points(elast_prop, src.zq);

    %%% Constant values %%%
    % Number of layers (INCLUDE halfspace)
    Nlayer = size(model_prop, 1);

    % Layers to record output at the top
    irec = (model_prop(:, end) == 1);

    % Calculate Lame parameters
    % Unit: GPa = 1e9 Pa
    mu = model_prop(:, 1) .* model_prop(:, 3).^2;
    lambda = model_prop(:, 1) .* model_prop(:, 2).^2 - 2.*mu;

    % Layer thickness & depth of layer top
    % Unit: km
    h = model_prop(:, 4);  ztop = model_prop(:, 5);

    %%% Mesh grid %%%
    % Spatial axes [km]
    Nx = length(src.xh);  dx = src.xh(2)-src.xh(1);
    Ny = length(src.yh);  dy = src.yh(2)-src.yh(1);

    %%% FFT parameters %%%
    % Wavenumber samples [rad/km]
    kx = 2*pi .* [0:Nx/2 (-Nx/2+1):-1]' ./ (Nx*dx);
    ky = 2*pi .* [0:Ny/2 (-Ny/2+1):-1]  ./ (Ny*dy);

    % Radial wavenumber samples [rad/km]
    Nr = max(Nx, Ny) * 10;  dr = min(dx, dy) / sqrt(3);
    kr = 2*pi .* (1:Nr/2)' ./ (Nr*dr);

    %%% Overflow guard %%%
    % The solution is propagated upward from the halfspace, so it grows like
    % exp(k*z) and the propagator entries overflow once k*z exceeds ~709. Past
    % that point every output field silently becomes NaN, so refuse instead.
    % k_max is set by the grid spacing, hence the suggested remedies below.
    if max(kr) * ztop(Nlayer) > 700
        error('solve_ds:tooDeep', ...
            ['Model is too deep for this grid: k_max*depth = %.0f exceeds 700 ' ...
             '(k_max = %.1f rad/km from dx = %g km, model depth = %g km).\n' ...
             'Coarsen the grid, or make the layered model shallower -- a mode ' ...
             'of wavenumber k cannot sense structure deeper than about 1/k, so ' ...
             'truncating the model below %.3g km changes nothing physically.'], ...
            max(kr)*ztop(Nlayer), max(kr), min(dx,dy), ztop(Nlayer), 700/max(kr));
    end

    %%% Initial homogeneous solution %%%
    lambda0 = lambda(end);  mu0 = mu(end);

    ds1j = [1/(2*mu0) ./ kr,   1/(2*mu0) ./ kr, ...
            ones(size(kr)),    ones(size(kr))]';

    ds2j = [(lambda0+2*mu0) ./ kr.^2 ./(2*mu0*(lambda0+mu0)), ...
            -1 ./ kr.^2 ./ (2*(lambda0+mu0)), ...
            1 ./ kr,           zeros(size(kr))]';

    %%% Propagator method %%%
    % Initialize arrays
    Nk = length(kr);
    ds1_surf = zeros(Nk, 4, Nlayer);
    ds2_surf = zeros(Nk, 4, Nlayer);

    % Loop over layers (including halfspace with h(end) = 0), vectorized over
    % wavenumber: the closed-form propagator below replaces a per-wavenumber expm
    for i = Nlayer:-1:1

        % Propagator matrix (closed form of expm(A*h), see psv_propagator)
        Pk = psv_propagator(kr, h(i), lambda(i), mu(i));

        ds1j = squeeze(pagemtimes(Pk, reshape(ds1j, 4, 1, Nk)));
        ds2j = squeeze(pagemtimes(Pk, reshape(ds2j, 4, 1, Nk)));

        % Record output
        ds1_surf(:,:,i) = ds1j';  ds2_surf(:,:,i) = ds2j';
    end


    %%% Output struct %%%
    output.xh = src.xh;  output.yh = src.yh;  output.dx = dx;  output.dy = dy;
    output.kx = kx;  output.ky = ky;  output.kr = kr;
    output.zq = model_prop(irec, 5);  output.prop = model_prop(irec, 1:3);
    output.ds1 = ds1_surf(:,:,irec);  output.ds2 = ds2_surf(:,:,irec);
end

%% Function: Closed-form propagator matrix for the P-SV system

% Exact expression for expm(A*h) with A as in Eq. (1) of the documentation,
% grouped in cosh/sinh so that nothing worse than cosh(k*h) can overflow.
% Returns a 4 x 4 x Nk array for the column vector of wavenumbers k. The layer
% thickness h may be a scalar or carry one value per wavenumber.

function Pk = psv_propagator(k, h, lambda, mu)

    k = k(:);  Nk = length(k);

    % Dimensionless modulus ratios
    sigma = lambda + 2*mu;
    al = mu/sigma;  be = (lambda+mu)/sigma;  ga = (lambda+3*mu)/sigma;

    % cosh & sinh of k*h
    x = k.*h;  ch = cosh(x);  sh = sinh(x);  km2 = 2*k*mu;

    Pk = zeros(4, 4, Nk);
    Pk(1,1,:) =  ch + be*x.*sh;            Pk(1,2,:) = -(al*sh + be*x.*ch);
    Pk(1,3,:) =  (ga*sh + be*x.*ch)./km2;  Pk(1,4,:) = -(be*x.*sh)./km2;
    Pk(2,1,:) =  be*x.*ch - al*sh;         Pk(2,2,:) =  ch - be*x.*sh;
    Pk(2,3,:) =  (be*x.*sh)./km2;          Pk(2,4,:) =  (ga*sh - be*x.*ch)./km2;
    Pk(3,1,:) =  be*km2.*(sh + x.*ch);     Pk(3,2,:) = -be*km2.*(x.*sh);
    Pk(3,3,:) =  ch + be*x.*sh;            Pk(3,4,:) =  al*sh - be*x.*ch;
    Pk(4,1,:) =  be*km2.*(x.*sh);          Pk(4,2,:) =  be*km2.*(sh - x.*ch);
    Pk(4,3,:) =  al*sh + be*x.*ch;         Pk(4,4,:) =  ch - be*x.*sh;
end

%% Function: Insert depth query points into layered model

function model_prop = insert_query_points(elast_prop, zq)

    %%% Ensure the last row has zero thickness (halfspace) %%%
    elast_prop(end, 4) = 0;
    
    %%% Combine depth sample points %%%
    % Depth of interfaces [km]
    z_layer = cumsum(elast_prop(1:end-1, 4));

    % Insert depth query points
    z_new = sort(union(zq(zq<z_layer(end) & zq>0), z_layer));


    %%% Initialize new layered model %%%
    model_prop = zeros(length(z_new)+1, 6);

    % 4th column: Layer thickness
    model_prop(2:end-1, 4) = diff(z_new);  model_prop(1, 4) = z_new(1);

    % 5th column: Depth of layer top
    model_prop(2:end, 5) = cumsum(model_prop(1:end-1, 4));

    % 6th (last) column: If to record output at the top
    model_prop(2:end, 6) = ismember(z_new, zq);  model_prop(1, 6) = 1;
    
    % Bottom halfspace
    model_prop(end, 1:4) = elast_prop(end, 1:4);

    % Assign each layer properties
    ind = discretize(z_new, union(0,z_layer), 'IncludedEdge', 'right');
    model_prop(1:end-1, 1:3) = elast_prop(ind, 1:3);

end
