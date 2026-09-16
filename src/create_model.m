%% Function: Create layered model from a continuous depth profile

% Each layer takes the profile value at its MIDPOINT, which makes the
% discretization second-order accurate in the layer thickness: halving the
% thickness cuts the error by four. Sampling at the layer top is only
% first-order, and for the same number of layers is far less accurate.

% This is for CONTINUOUS profiles. A model with genuine interfaces is better
% written as elast_prop directly, so that a node sits exactly on each boundary:
% that is exact and costs nothing, whereas discretizing a step here leaves the
% interface misplaced by up to one layer thickness unless it happens to land on
% a node. For a mixed model, call this once per gradient zone and stack the
% blocks, so every real interface still gets its own node.

% Input profile is either a struct of function handles (rho, vp, vs) of
% depth [km], or a matrix [depth, rho, vp, vs] to be interpolated.
% H is the depth of the base of the profile [km], below which the halfspace
% lies. Nlayer is the number of finite layers, so the output has Nlayer+1 rows.

% Optional param fields:
%   grading   Layer interfaces at H*(j/Nlayer)^grading. This redistributes a
%             fixed number of layers, it does not change how many there are.
%             1 gives uniform layers, larger values thin them toward the
%             surface. Worth using when the profile is steepest near the
%             surface, as shallow velocity profiles usually are: for
%             vs ~ z^0.3 a uniform model converges only as Nlayer^-1.3, while
%             grading = 3 restores second order and is 20 to 80 times more
%             accurate at the same layer count. For a near-linear profile it
%             makes little difference, so the default is 1.
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
