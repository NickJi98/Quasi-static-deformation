%% Function: Quasi-static seismic modeling in layered medium (P-SV system)

function [sol_pm, ds_pm] = qs_model(src, elast_prop)

    % Depth query points
    if ~src.depth_query
        src.zq = 0;
    end

    % Default zero surface loads for missing fields
    % (pp: normal pressure; tx, ty: horizontal traction along x, y)
    if ~isfield(src, 'pp');  src.pp = zeros(length(src.xh), length(src.yh));  end
    if ~isfield(src, 'tx');  src.tx = zeros(length(src.xh), length(src.yh));  end
    if ~isfield(src, 'ty');  src.ty = zeros(length(src.xh), length(src.yh));  end

    % Solve displacement-stress vector toward the surface
    ds_pm = solve_ds(src, elast_prop);

    % For static loading, add the third axis
    if ismatrix(src.pp)
        src.pp = src.pp(:,:,1);
    end
    
    % Numerical solution
    sol_pm = calc_layer(src, ds_pm);

    % Return modeling results
    % (For vertical displacement, positive downward)
    % (For stress component, positive for tensile direction)
    sol_pm.uz = -sol_pm.uz;

    % Add time info to output
    % (a static load is just a time axis with one sample, so the time vector is
    %  carried through whenever it is supplied, independent of the load length)
    if isfield(src, 'time')
        sol_pm.time = src.time;
    end
end
