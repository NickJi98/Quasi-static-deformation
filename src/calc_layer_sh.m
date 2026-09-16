%% Function: Calculate displacement & stress for layered medium (SH system)

function output = calc_layer_sh(src, sol_ds)

    %%% Read input struct %%%
    Nx = length(sol_ds.xh);  Ny = length(sol_ds.yh);  Nz = length(sol_ds.zq);
    dx = sol_ds.dx;  dy = sol_ds.dy;
    kx = sol_ds.kx;  ky = sol_ds.ky;  kr = sol_ds.kr;
    dsh = sol_ds.dsh;
    stress_flag = src.include_stress;


    %%% Convert to Cartesian coordinate %%%
    % Wavenumber grid [rad/km]
    [Kx, Ky] = meshgrid(kx, ky);  Kx = Kx';  Ky = Ky';
    Kr = sqrt(Kx.^2 + Ky.^2);     clear Kx Ky;

    % Initialize arrays
    dsh_xy = zeros(Nx,Ny,Nz,2);

    % Interpolation
    for j = 1:2
        for iz = 1:Nz
            dsh_xy(:,:,iz,j) = interp1(kr, dsh(:,j,iz), Kr, 'linear', 0);
        end
    end

    % Fix k = 0 component
    Kr(1, 1) = Inf;


    %%% Solve linear system at each time %%%
    % Number of time steps
    Nt = max(size(src.tx, 3), size(src.ty, 3));

    % Horizontal traction (Positive along +x, +y direction)
    fk_tx = fft(fft(src.tx,[],1),[],2) .* dx*dy;
    fk_ty = fft(fft(src.ty,[],1),[],2) .* dx*dy;

    % Surface value of V2 from transverse traction (-k*V2 = i*kx*syz - i*ky*sxz)
    fk_V2s = -1j .* (kx.*fk_ty - ky.*fk_tx) ./ Kr;

    % Initialize output arrays
    uz = zeros(Nx,Ny,Nz,Nt);  ux = zeros(Nx,Ny,Nz,Nt);  uy = zeros(Nx,Ny,Nz,Nt);
    if stress_flag
        sxx = zeros(Nx,Ny,Nz,Nt);  syy = zeros(Nx,Ny,Nz,Nt);  szz = zeros(Nx,Ny,Nz,Nt);
        sxy = zeros(Nx,Ny,Nz,Nt);  sxz = zeros(Nx,Ny,Nz,Nt);  syz = zeros(Nx,Ny,Nz,Nt);
    end

    % Elastic properties for calculating stress
    prop = sol_ds.prop;
    mu = reshape(prop(:, 1).*prop(:, 3).^2, [1,1,Nz]);

    % Reciprocal of the (time-independent) SH denominator, precomputed once:
    % V2 = c3 * dsh(2) at the surface, so c3 = V2 / dsh(2). Only the right-hand
    % side changes with time, which matters for quasi-static runs with many steps.
    ivh = 1 ./ dsh_xy(:,:,1,2);

    % Drop wavenumbers with no solution: k = 0 (REMOVED, as before) and any
    % singular entry, so they contribute zero instead of spreading NaN over the
    % whole field through the inverse FFT
    ivh(~isfinite(ivh)) = 0;  ivh(1,1) = 0;

    % Time samples of the traction load
    Nt_ts = size(fk_V2s, 3);

    for it = 1:Nt
        % Coefficient of the SH system (single homogeneous solution)
        c = fk_V2s(:,:,min(it,Nt_ts)) .* ivh;

        % Surface displacement
        fk_V1 = c .* dsh_xy(:,:,:,1);

        % Solve horizontal displacement (uz = 0 for SH system)
        fk_ux = -1j.*fk_V1.*ky./Kr;  fk_uy = 1j.*fk_V1.*kx./Kr;

        ux(:,:,:,it) = real(ifft(ifft(fk_ux,[],1),[],2)) ./ (dx*dy);
        uy(:,:,:,it) = real(ifft(ifft(fk_uy,[],1),[],2)) ./ (dx*dy);

        if stress_flag
            % Surface stress
            fk_V2 = c .* dsh_xy(:,:,:,2);

            % Solve shear stress (szz = 0 for SH system)
            fk_sxz = -1j.*fk_V2.*ky./Kr;  fk_syz = 1j.*fk_V2.*kx./Kr;
            fk_sxy = mu .* 1j.*(ky.*fk_ux + kx.*fk_uy);

            % Solve normal stress (zero dilatation for SH system)
            fk_sxx = 2.*mu .* 1j.*kx.*fk_ux;
            fk_syy = 2.*mu .* 1j.*ky.*fk_uy;

            sxx(:,:,:,it) = real(ifft(ifft(fk_sxx,[],1),[],2)) ./ (dx*dy);
            syy(:,:,:,it) = real(ifft(ifft(fk_syy,[],1),[],2)) ./ (dx*dy);
            sxy(:,:,:,it) = real(ifft(ifft(fk_sxy,[],1),[],2)) ./ (dx*dy);
            sxz(:,:,:,it) = real(ifft(ifft(fk_sxz,[],1),[],2)) ./ (dx*dy);
            syz(:,:,:,it) = real(ifft(ifft(fk_syz,[],1),[],2)) ./ (dx*dy);
        end
    end


    %%% Output struct %%%
    output.x = sol_ds.xh;  output.y = sol_ds.yh;  output.zq = sol_ds.zq;
    output.ux = ux;  output.uy = uy;  output.uz = uz;

    if stress_flag
        output.sxx = sxx;  output.syy = syy;  output.szz = szz;
        output.sxy = sxy;  output.sxz = sxz;  output.syz = syz;
    end

end
