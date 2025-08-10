%% Function: Read CM1 output data

function data_struct = create_simple_src(cm1_file, Uc)

    %%% Create data structure

    % Mesh grid
    data_struct.xh = double(ncread(cm1_file, 'xh'));
    data_struct.yh = double(ncread(cm1_file, 'yh'));
    data_struct.time = ncread(cm1_file, 'time');

    % FFT wavenumber samples
    xh = data_struct.xh;  yh = data_struct.yh;  t = 0:0.5:1800-0.5;
    Nx = length(xh);  dx = xh(2) - xh(1);
    k = (1:Nx/2) .* (2*pi/(Nx*dx));  k = reshape(k, [1,1,1,Nx/2]);

    % Amplitude spectrum
    w = Uc .* k;  C = 3.5;  tau = 2;
    % amp = 0.811 * C^2 / tau .* (1 + 0.1792 .* (w .* tau).^2).^(-7/6);
    % loglog(squeeze(2*pi./w), squeeze(amp), 'k-');

    % Random phase
    phi = 2 * pi * rand(1, Nx/2);  phi = reshape(phi, [1,1,1,Nx/2]);

    % Surface pressure [Pa]
    xh = reshape(xh, [Nx, 1, 1]);  yh = reshape(yh, [1, Nx, 1]);
    t = reshape(t, [1, 1, length(t)]);
    pp = sum(cos(k.*xh - w.*t + phi), 4);
    pp = pp ./ max(pp, [], 'all') .* 14.14 * 5;
    pp = repmat(pp, [1,Nx,1]);
    data_struct.pp = pp;
    
end
