%% Function: Sorrells-type solution of displacement (shear traction wave)

function output = calc_sorrells_shear(src, zq, elast_prop)

    %%% Elastic properties %%%
    rho = elast_prop(1);  vp = elast_prop(2);  vs = elast_prop(3);

    % Shear modulus [GPa]
    mu = rho * vs^2;

    %%% Traction wave %%%
    T0 = src.amp;               % Amplitude [Pa]
    kw = 2*pi/src.wavelen;      % Wavenumber [rad/km]
    az = src.waveaz;            % Wave azimuth [rad]
    taz = src.taz;              % Traction azimuth [rad]

    % Traction components along & transverse to the wave direction
    T_par = T0 * cos(taz - az);  T_perp = T0 * sin(taz - az);

    %%% Sorrells-type solution %%%
    % Vertical displacement: Positive for upward motion

    % P-SV response to in-plane traction [μm]
    % (upar in phase with traction, uz shifted by 90°)
    upar = T_par/(2*mu*kw) .* (vp^2/(vp^2-vs^2) - kw.*zq) .* exp(-kw.*zq);
    uz   = T_par/(2*mu*kw) .* (vs^2/(vp^2-vs^2) + kw.*zq) .* exp(-kw.*zq);

    % SH response to transverse traction [μm]
    % (uperp in phase with traction)
    uperp = T_perp/(mu*kw) .* exp(-kw.*zq);

    % Cartesian components (upar & uperp share the same phase)
    ux = upar .* cos(az) - uperp .* sin(az);
    uy = upar .* sin(az) + uperp .* cos(az);

    %%% Output struct %%%
    output.uz = abs(uz);     output.ux = abs(ux);     output.uy = abs(uy);
    output.zq = zq;

end
