%% Function: Sorrells solution of displacement

function output = calc_sorrells(src, zq, elast_prop)

    %%% Elastic properties %%%
    rho = elast_prop(1);  vp = elast_prop(2);  vs = elast_prop(3);

    % Shear modulus [GPa]
    mu = rho * vs^2;

    %%% Pressure wave %%%
    P0 = src.amp;               % Amplitude [Pa]
    kw = 2*pi/src.wavelen;      % Wavenumber [rad/km]
    az = src.waveaz;            % Azimuth [rad]

    %%% Sorrells solution %%%
    % Vertical displacement: Positive for upward motion

    % Displacement amplitude [μm]
    uz = -P0/(2*mu*kw) .* (vp^2/(vp^2-vs^2) + kw.*zq) .* exp(-kw.*zq);
    uh =  P0/(2*mu*kw) .* (vs^2/(vp^2-vs^2) - kw.*zq) .* exp(-kw.*zq);
    ux = uh * cos(az);  uy = uh * sin(az);

    %%% Output struct %%%
    output.uz = abs(uz);     output.ux = abs(ux);     output.uy = abs(uy);
    output.zq = zq;

end
