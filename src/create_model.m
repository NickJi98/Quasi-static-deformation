%% Function: Discretize a continuous subsurface profile into layers

% Build the elast_prop matrix that qs_model expects from a continuously
% varying profile. Each layer takes the profile value at its MIDPOINT, which
% makes the discretization second-order accurate in the layer thickness:
% halving the thickness cuts the error by four. Sampling at the layer top
% instead is only first-order, so for the same number of layers the midpoint
% rule is typically two to three orders of magnitude more accurate.

% Usage
%   elast_prop = create_model(profile, H, Nlayer)
%   elast_prop = create_model(profile, H, Nlayer, 'grading', 2)
%
%   profile : struct with fields rho, vp, vs, each a function handle of
%             depth [km]  (e.g. profile.vs = @(z) 0.27 + 0.63*z)
%             OR a matrix [depth, rho, vp, vs] that is interpolated
%   H       : depth of the base of the profile [km]
%   Nlayer  : number of finite layers (the output has Nlayer+1 rows, the
%             last being the underlying halfspace)
%
% Options
%   'grading'   g > 0. Layer interfaces are placed at H*(j/Nlayer)^g, so
%               g = 1 gives uniform layers and g > 1 makes them thinner near
%               the surface, where short-wavelength loads concentrate.
%   'halfspace' [rho vp vs] for the underlying halfspace. Defaults to the
%               profile evaluated at depth H.
%   'method'    interpolation method for a tabulated profile ('pchip').

function elast_prop = create_model(profile, H, Nlayer, varargin)

    %%% Parse options %%%
    opt = struct('grading', 1, 'halfspace', [], 'method', 'pchip');
    for i = 1:2:numel(varargin)
        name = lower(varargin{i});
        if ~isfield(opt, name)
            error('create_model:badOption', 'Unknown option ''%s''.', varargin{i});
        end
        opt.(name) = varargin{i+1};
    end

    if ~isscalar(H) || H <= 0
        error('create_model:badDepth', 'H must be a positive scalar depth [km].');
    end
    if ~isscalar(Nlayer) || Nlayer < 1 || mod(Nlayer,1) ~= 0
        error('create_model:badNlayer', 'Nlayer must be a positive integer.');
    end
    if ~isscalar(opt.grading) || opt.grading <= 0
        error('create_model:badGrading', 'grading must be a positive scalar.');
    end

    %%% Profile evaluator %%%
    fprop = make_evaluator(profile, H, opt.method);

    %%% Layer interfaces & midpoints [km] %%%
    % Graded spacing: g = 1 uniform, g > 1 refines toward the surface
    zface = H .* ((0:Nlayer)' ./ Nlayer).^opt.grading;
    hlayer = diff(zface);
    zmid = zface(1:end-1) + hlayer./2;

    %%% Assign layer properties at the midpoints %%%
    elast_prop = zeros(Nlayer+1, 4);
    elast_prop(1:Nlayer, 1:3) = fprop(zmid);
    elast_prop(1:Nlayer, 4)   = hlayer;

    %%% Underlying halfspace (zero thickness) %%%
    if isempty(opt.halfspace)
        elast_prop(end, 1:3) = fprop(H);
    else
        hs = opt.halfspace(:)';
        if numel(hs) ~= 3
            error('create_model:badHalfspace', 'halfspace must be [rho vp vs].');
        end
        elast_prop(end, 1:3) = hs;
    end
    elast_prop(end, 4) = 0;

    %%% Sanity check %%%
    if any(~isfinite(elast_prop(:))) || any(elast_prop(:,1:3) <= 0, 'all')
        error('create_model:badProfile', ...
            'Profile produced non-positive or non-finite rho, vp or vs.');
    end
    if any(elast_prop(1:Nlayer,2) <= elast_prop(1:Nlayer,3).*sqrt(2), 'all')
        warning('create_model:poisson', ...
            'vp <= sqrt(2)*vs in some layers, i.e. negative Poisson ratio.');
    end
end

%% Function: Build a callable [rho vp vs] evaluator from the input profile

function fprop = make_evaluator(profile, H, method)

    if isstruct(profile)
        need = {'rho', 'vp', 'vs'};
        if ~all(isfield(profile, need))
            error('create_model:badProfile', ...
                'Profile struct needs fields rho, vp and vs.');
        end
        if ~all(cellfun(@(f) isa(profile.(f), 'function_handle'), need))
            error('create_model:badProfile', ...
                'Profile struct fields must be function handles of depth.');
        end
        fprop = @(z) [profile.rho(z(:)).*ones(numel(z),1), ...
                      profile.vp(z(:)).*ones(numel(z),1), ...
                      profile.vs(z(:)).*ones(numel(z),1)];

    elseif isnumeric(profile) && size(profile,2) == 4
        zt = profile(:,1);
        if min(zt) > 0 || max(zt) < H
            error('create_model:badProfile', ...
                ['Tabulated profile spans %g to %g km but must cover 0 to %g km ' ...
                 '(no extrapolation).'], min(zt), max(zt), H);
        end
        fprop = @(z) interp1(zt, profile(:,2:4), z(:), method);

    else
        error('create_model:badProfile', ...
            ['Profile must be a struct of function handles (rho, vp, vs) ' ...
             'or a matrix [depth, rho, vp, vs].']);
    end
end
