%% Function: Generate surface horizontal traction field

function data_struct = create_tsfc(test_case, mesh, param)

    %%% Generate mesh grid %%%
    xh = (0:mesh.Nx-1)' .* mesh.dx;
    try
        yh = (0:mesh.Ny-1) .* mesh.dy;
    catch
        yh = xh';   % If y not provided, consider square domain
    end
    data_struct.xh = xh;  data_struct.yh = yh;

    %%% Create surface traction magnitude %%%
    switch test_case

        % Gaussian load
        case 1
            tt = gaussian_tsfc(xh, yh, param);
            data_struct.sigma = sqrt(param.sigma_x * param.sigma_y);

        % Wave snapshot
        case 2
            [tt, time, wavelen, wavespd, waveaz] = period_tsfc(xh, yh, mesh, param);
            data_struct.time = time;  data_struct.amp = param.amp;
            data_struct.wavelen = wavelen;  data_struct.wavespd = wavespd;
            data_struct.waveaz = waveaz;

        % White spectrum
        case 4
            tt = white_tsfc(mesh);

    end

    %%% Traction components along the given azimuth %%%
    % (taz: traction azimuth [rad] measured from +x axis)
    data_struct.taz = param.taz;
    data_struct.tx = tt .* cos(param.taz);
    data_struct.ty = tt .* sin(param.taz);

end

%% Function: Gaussian load

function tt = gaussian_tsfc(xh, yh, param)

    % Parameters (amp: Pa, sigma: km)
    amp = param.amp;  sigma_x = param.sigma_x;  sigma_y = param.sigma_y;
    x0 = param.x0;    y0 = param.y0;

    % Gaussian source
    tt = amp .* exp(-(xh-x0).^2./(2*sigma_x^2) - (yh-y0).^2./(2*sigma_y^2)) ...
      ./ (2*pi*sigma_x*sigma_y);

end

%% Function: Periodic load

function [tt, time, wavelen, wavespd, waveaz] = period_tsfc(xh, yh, mesh, param)

    % Mesh wavenumber (kx, ky: rad/km)
    kx = 2*pi / (mesh.Nx * mesh.dx);
    try
        ky = 2*pi / (mesh.Ny * mesh.dy);
    catch
        ky = kx;
    end

    % Time mesh (t: s, two cycles)
    time = linspace(0, 2/param.fw, 16);

    % Parameters (amp: Pa, kw: rad/km, omega: rad/s)
    amp = param.amp;
    kw_x = param.Nw_x * kx;  kw_y = param.Nw_y * ky;
    kw = sqrt(kw_x^2 + kw_y^2);
    omega = 2*pi * param.fw;

    % Traction wave
    tt = amp .* cos(kw_x.*xh + kw_y.*yh - ...
        omega.*reshape(time, [1,1,length(time)]));

    % Wavelength (km) & wave speed (km/s)
    wavelen = 2*pi/kw;  wavespd = wavelen * param.fw;
    waveaz = atan(kw_y/kw_x);

end

%% Function: White spectrum load

function tt = white_tsfc(mesh)

    % Mesh wavenumber (kx, ky: rad/km)
    Nx = mesh.Nx;  dx = mesh.dx;
    kx = 2*pi .* [0:Nx/2 (-Nx/2+1):-1]' ./ (Nx*dx);
    try
        Ny = mesh.Ny;  dy = mesh.dy;
        ky = 2*pi .* [0:Ny/2 (-Ny/2+1):-1] ./ (Ny*dy);
    catch
        ky = kx';  Ny = Nx;  dy = dx;
    end

    % White spectrum
    fk_tt = ones(Nx, Ny);

    % Shift to domain center
    fk_tt = fk_tt .* exp(-1j.*kx .* Nx*dx/2) .* exp(-1j.*ky .* Ny*dy/2);

    % Inverse FFT
    tt = real(ifft(ifft(fk_tt,[],1),[],2)) ./ (dx*dy);

end
