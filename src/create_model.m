%% Function: Create layered model from a continuous depth profile

% Each layer takes the profile value at its MIDPOINT, which makes the
% discretization second-order accurate in the layer thickness: halving the
% thickness cuts the error by four. Sampling at the layer top is only
% first-order, and for the same number of layers is far less accurate.

% Input profile is either a struct of function handles (rho, vp, vs) of
% depth [km], or a matrix [depth, rho, vp, vs] to be interpolated.
% H is the depth of the base of the profile [km], below which the halfspace
% lies. Nlayer is the number of finite layers, so the output has Nlayer+1 rows.

% Optional param fields:
%   grading   Layer interfaces at H*(j/Nlayer)^grading, so 1 gives uniform
%             layers and larger values thin them toward the surface
%   halfspace [rho vp vs] of the bottom halfspace, default is profile at H
%   method    Interpolation method for a tabulated profile, default 'pchip'

function elast_prop = create_model(profile, H, Nlayer, param)

    %%% Optional parameters %%%
    try
        grading = param.grading;
    catch
        grading = 1;    % Uniform layer thickness
    end
    try
        method = param.method;
    catch
        method = 'pchip';
    end

    %%% Layer interfaces & midpoints %%%
    % Unit: km
    zface = H .* ((0:Nlayer)' ./ Nlayer).^grading;
    hlayer = diff(zface);
    zmid = zface(1:end-1) + hlayer./2;

    %%% Assign layer properties at the midpoints %%%
    % Units: g/cm^3, km/s, km/s, km
    elast_prop = zeros(Nlayer+1, 4);
    elast_prop(1:Nlayer, 1:3) = eval_profile(profile, zmid, method);
    elast_prop(1:Nlayer, 4) = hlayer;

    %%% Bottom halfspace (zero thickness) %%%
    try
        elast_prop(end, 1:3) = param.halfspace;
    catch
        elast_prop(end, 1:3) = eval_profile(profile, H, method);
    end

    %%% Check the resulting model %%%
    if any(~isfinite(elast_prop), 'all') || any(elast_prop(:,1:3) <= 0, 'all')
        error('Profile gives non-positive or non-finite rho, vp or vs!');
    end

end

%% Function: Evaluate the profile at given depths

function prop = eval_profile(profile, zq, method)

    zq = zq(:);

    % Struct of function handles
    if isstruct(profile)
        prop = [profile.rho(zq).*ones(size(zq)), ...
                profile.vp(zq) .*ones(size(zq)), ...
                profile.vs(zq) .*ones(size(zq))];

    % Tabulated profile [depth, rho, vp, vs]
    elseif isnumeric(profile) && size(profile, 2) == 4
        prop = interp1(profile(:,1), profile(:,2:4), zq, method);

    else
        error('Profile should be a struct of rho, vp, vs or a 4-column matrix!');
    end

end
