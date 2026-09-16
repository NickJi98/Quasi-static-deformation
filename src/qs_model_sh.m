%% Function: Quasi-static seismic modeling in layered medium (SH system)

function [sol_pm, ds_pm] = qs_model_sh(src, elast_prop)

    % Depth query points
    if ~src.depth_query
        src.zq = 0;
    end

    % Default zero surface loads for missing fields
    % (tx, ty: horizontal traction along x, y)
    if ~isfield(src, 'tx');  src.tx = zeros(length(src.xh), length(src.yh));  end
    if ~isfield(src, 'ty');  src.ty = zeros(length(src.xh), length(src.yh));  end

    % Solve displacement-stress vector toward the surface
    ds_pm = solve_ds_sh(src, elast_prop);

    % For static loading, add the third axis
    if ismatrix(src.tx);  src.tx = src.tx(:,:,1);  end
    if ismatrix(src.ty);  src.ty = src.ty(:,:,1);  end

    % Numerical solution
    sol_pm = calc_layer_sh(src, ds_pm);

    % Return modeling results
    % (SH system has no vertical displacement: uz = 0, szz = 0)
    % (For stress component, positive for tensile direction)

    % Add time info to output
    % (a static load is just a time axis with one sample, so the time vector is
    %  carried through whenever it is supplied, independent of the load length)
    if isfield(src, 'time')
        sol_pm.time = src.time;
    end
end
