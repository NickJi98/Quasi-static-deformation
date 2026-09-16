%% Function: Solve displacement-stress vector (SH system)

% Return output vectors at different depths specified by zq
% Output at surface z = 0 is always recorded

function output = solve_ds_sh(src, elast_prop)

    %%% Insert depth query points into layered model %%%
    model_prop = insert_query_points(elast_prop, src.zq);

    %%% Constant values %%%
    % Number of layers (INCLUDE halfspace)
    Nlayer = size(model_prop, 1);

    % Layers to record output at the top
    irec = (model_prop(:, end) == 1);

    % Calculate shear modulus
    % Unit: GPa = 1e9 Pa
    mu = model_prop(:, 1) .* model_prop(:, 3).^2;

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
    % See solve_ds.m: the upward-propagated solution grows like exp(k*z) and
    % overflows past k*z ~ 709, which would silently turn every field into NaN.
    if max(kr) * ztop(Nlayer) > 700
        error(['Model is too deep for this grid: k_max*depth = %.0f exceeds 700! ' ...
            '(k_max = %.1f rad/km from dx = %g km, model depth = %g km) ' ...
            'Coarsen the grid, or make the layered model shallower. A mode of ' ...
            'wavenumber k cannot sense structure deeper than about 1/k, so ' ...
            'truncating the model below %.3g km changes nothing physically.'], ...
            max(kr)*ztop(Nlayer), max(kr), min(dx,dy), ztop(Nlayer), 700/max(kr));
    end

    %%% Initial homogeneous solution %%%
    mu0 = mu(end);

    dshj = [1/mu0 ./ kr,   ones(size(kr))]';

    %%% Propagator method %%%
    % Initialize arrays
    Nk = length(kr);
    dsh_surf = zeros(Nk, 2, Nlayer);

    % Loop over layers (including halfspace with h(end) = 0), vectorized over
    % wavenumber: the closed-form propagator below replaces a per-wavenumber expm
    for i = Nlayer:-1:1

        % Propagator matrix (closed form of expm(B*h))
        Pk = sh_propagator(kr, h(i), mu(i));

        dshj = squeeze(pagemtimes(Pk, reshape(dshj, 2, 1, Nk)));

        % Record output
        dsh_surf(:,:,i) = dshj';
    end


    %%% Output struct %%%
    output.xh = src.xh;  output.yh = src.yh;  output.dx = dx;  output.dy = dy;
    output.kx = kx;  output.ky = ky;  output.kr = kr;
    output.zq = model_prop(irec, 5);  output.prop = model_prop(irec, 1:3);
    output.dsh = dsh_surf(:,:,irec);
end

%% Function: Closed-form propagator matrix for the SH system

% Exact expression for expm(B*h) with B as in the SH section of the
% documentation. Returns a 2 x 2 x Nk array for the wavenumber vector k; the
% thickness h may be a scalar or carry one value per wavenumber.

function Pk = sh_propagator(k, h, mu)

    k = k(:);  Nk = length(k);

    % cosh & sinh of k*h
    x = k.*h;  ch = cosh(x);  sh = sinh(x);

    Pk = zeros(2, 2, Nk);
    Pk(1,1,:) = ch;               Pk(1,2,:) = sh ./ (mu.*k);
    Pk(2,1,:) = mu.*k .* sh;      Pk(2,2,:) = ch;
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
