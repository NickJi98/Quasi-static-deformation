%% Function: Generate surface pressure field

function data_struct = create_psfc(test_case, mesh, param)

    %%% Generate mesh grid %%%
    xh = (0:mesh.Nx-1)' .* mesh.dx;
    try
        yh = (0:mesh.Ny-1) .* mesh.dy;
    catch
        yh = xh';   % If y not provided, consider square domain
    end
    data_struct.xh = xh;  data_struct.yh = yh;

    %%% Create surface pressure %%%
    switch test_case

        % Gaussian load
        case 1
            pp = gaussian_psfc(xh, yh, param);
            data_struct.sigma = sqrt(param.sigma_x * param.sigma_y);
        
        % Wave snapshot
        case 2
            [pp, time, wavelen, wavespd, waveaz] = period_psfc(xh, yh, mesh, param);
            data_struct.time = time;  data_struct.amp = param.amp;
            data_struct.wavelen = wavelen;  data_struct.wavespd = wavespd;
            data_struct.waveaz = waveaz;

        % CM1 snapshot
        case 3
            [pp, xh, time] = cm1_psfc();
            data_struct.xh = xh;  data_struct.yh = xh';
            data_struct.time = time;

        % White spectrum
        case 4
            pp = white_psfc(mesh);
        
    end

    %%% NOT remove the mean of pressure field %%%
    data_struct.pp = pp;

end

%% Function: Gaussian load

function pp = gaussian_psfc(xh, yh, param)

    % Parameters (amp: Pa, sigma: km)
    amp = param.amp;  sigma_x = param.sigma_x;  sigma_y = param.sigma_y;
    x0 = param.x0;    y0 = param.y0;

    % Gaussian source
    pp = amp .* exp(-(xh-x0).^2./(2*sigma_x^2) - (yh-y0).^2./(2*sigma_y^2)) ...
      ./ (2*pi*sigma_x*sigma_y);

end

%% Function: Periodic load

function [pp, time, wavelen, wavespd, waveaz] = period_psfc(xh, yh, mesh, param)

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

    % Gaussian source
    pp = amp .* cos(kw_x.*xh + kw_y.*yh - ...
        omega.*reshape(time, [1,1,length(time)]));

    % Wavelength (km) & wave speed (km/s)
    wavelen = 2*pi/kw;  wavespd = wavelen * param.fw;
    waveaz = atan(kw_y/kw_x);

end

%% Function: Example CM1 surface pressure

function [pp, xh, time] = cm1_psfc()

    % CM1 output file
    src_dir = fileparts(mfilename('fullpath'));
    cm1_file = fullfile(src_dir, '..', 'data', 'cm1out.nc');

    % Read CM1 output (30 frames)
    xh = double(ncread(cm1_file, 'xh'));
    pp = double(squeeze(ncread(cm1_file, 'psfc', [1,1,900], [Inf,Inf,30])));
    time = ncread(cm1_file, 'time', 900, 30);

    % Remove spatial mean at each time
    pp = pp - mean(pp, [1,2]);

end

%% Function: White spectrum load

function pp = white_psfc(mesh)

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
    fk_pp = ones(Nx, Ny);

    % Shift to domain center
    fk_pp = fk_pp .* exp(-1j.*kx .* Nx*dx/2) .* exp(-1j.*ky .* Ny*dy/2);

    % Inverse FFT
    pp = real(ifft(ifft(fk_pp,[],1),[],2)) ./ (dx*dy);

end
