%% Function: Calculate displacement & stress for layered medium (P-SV system)

function output = calc_layer(src, sol_ds)

    %%% Read input struct %%%
    Nx = length(sol_ds.xh);  Ny = length(sol_ds.yh);  Nz = length(sol_ds.zq);
    dx = sol_ds.dx;  dy = sol_ds.dy;
    kx = sol_ds.kx;  ky = sol_ds.ky;  kr = sol_ds.kr;
    ds1 = sol_ds.ds1;   ds2 = sol_ds.ds2;
    stress_flag = src.include_stress;


    %%% Convert to Cartesian coordinate %%%
    % Wavenumber grid [rad/km]
    [Kx, Ky] = meshgrid(kx, ky);  Kx = Kx';  Ky = Ky';
    Kr = sqrt(Kx.^2 + Ky.^2);     clear Kx Ky;
    
    % Initialize arrays
    ds1_xy = zeros(Nx,Ny,Nz,4);  ds2_xy = zeros(Nx,Ny,Nz,4);
    
    % Interpolation
    parfor j = 1:4
        for iz = 1:Nz
            ds1_xy(:,:,iz,j) = interp1(kr, ds1(:,j,iz), Kr, 'linear', 0);
            ds2_xy(:,:,iz,j) = interp1(kr, ds2(:,j,iz), Kr, 'linear', 0);
        end
    end

    % Fix k = 0 component
    Kr(1, 1) = Inf;


    %%% Solve linear system at each time %%%
    % Number of time steps
    Nt = size(src.pp, 3);

    % Normal traction (Positive for tensile)
    fk_pp = -fft(fft(src.pp,[],1),[],2) .* dx*dy;

    % Matrix of the linear system (pagewise)
    A2 = cat(4, squeeze(ds1_xy(:,:,1,3:4)), squeeze(ds2_xy(:,:,1,3:4)));

    % Initialize output arrays
    uz = zeros(Nx,Ny,Nz,Nt);  ux = zeros(Nx,Ny,Nz,Nt);  uy = zeros(Nx,Ny,Nz,Nt);
    if stress_flag
        sxx = zeros(Nx,Ny,Nz,Nt);  syy = zeros(Nx,Ny,Nz,Nt);  szz = zeros(Nx,Ny,Nz,Nt);
        sxy = zeros(Nx,Ny,Nz,Nt);  sxz = zeros(Nx,Ny,Nz,Nt);  syz = zeros(Nx,Ny,Nz,Nt);
        
        % Elastic properties for calculating stress
        prop = sol_ds.prop;
        lambda = reshape(prop(:, 1).*prop(:, 2).^2, [1,1,Nz]);
        mu = reshape(prop(:, 1).*prop(:, 3).^2, [1,1,Nz]);
    end

    parfor it = 1:Nt
        % Vector of the linear system (pagewise)
        b2 = zeros(1, Nx, Ny, 2);
        b2(1,:,:,:) = cat(3, zeros(Nx, Ny), fk_pp(:,:,it));

        % Solve linear system (pagewise)
        c = pagemldivide(permute(A2, [3,4,1,2]), permute(b2, [4,1,2,3]));
        c = squeeze(permute(c, [3,4,1,2]));

        % Remove NaN (REMOVE k = 0 component)
        c(1,1,:) = 0;

        % Surface displacement
        fk_uz = c(:,:,1) .* ds1_xy(:,:,:,2) + c(:,:,2) .* ds2_xy(:,:,:,2);
        fk_U1 = c(:,:,1) .* ds1_xy(:,:,:,1) + c(:,:,2) .* ds2_xy(:,:,:,1);

        % Solve horizontal displacement
        fk_ux = 1j.*fk_U1.*kx./Kr;  fk_uy = 1j.*fk_U1.*ky./Kr;

        ux(:,:,:,it) = real(ifft(ifft(fk_ux,[],1),[],2)) ./ (dx*dy);
        uy(:,:,:,it) = real(ifft(ifft(fk_uy,[],1),[],2)) ./ (dx*dy);
        uz(:,:,:,it) = real(ifft(ifft(fk_uz,[],1),[],2)) ./ (dx*dy);

        if stress_flag
            % Surface stress
            fk_szz = c(:,:,1) .* ds1_xy(:,:,:,4) + c(:,:,2) .* ds2_xy(:,:,:,4);
            fk_U3  = c(:,:,1) .* ds1_xy(:,:,:,3) + c(:,:,2) .* ds2_xy(:,:,:,3);

            % Solve shear stress
            fk_sxz = 1j.*fk_U3.*kx./Kr;  fk_syz = 1j.*fk_U3.*ky./Kr;
            fk_sxy = mu .* 1j.*(ky.*fk_ux + kx.*fk_uy);
    
            % Solve normal stress
            fk_sxx = (lambda.*fk_szz + 1j.*(4.*mu.*(lambda+mu).*kx.*fk_ux ...
                + 2.*lambda.*mu.*ky.*fk_uy)) ./ (lambda+2.*mu);
            fk_syy = (lambda.*fk_szz + 1j.*(4.*mu.*(lambda+mu).*ky.*fk_uy ...
                + 2.*lambda.*mu.*kx.*fk_ux)) ./ (lambda+2.*mu);

            sxx(:,:,:,it) = real(ifft(ifft(fk_sxx,[],1),[],2)) ./ (dx*dy);
            syy(:,:,:,it) = real(ifft(ifft(fk_syy,[],1),[],2)) ./ (dx*dy);
            szz(:,:,:,it) = real(ifft(ifft(fk_szz,[],1),[],2)) ./ (dx*dy);
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
