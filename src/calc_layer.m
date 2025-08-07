%% Function: Calculate displacement & stress in layered medium (P-SV system)

function output = calc_layer(sol_ds, coeff, coord, include_stress)

    %%% Read input struct %%%
    Nx = length(sol_ds.xh);  Ny = length(sol_ds.yh);  Nz = length(sol_ds.zq);
    dx = sol_ds.dx;  dy = sol_ds.dy;
    kx = sol_ds.kx;  ky = sol_ds.ky;  kr = sol_ds.kr;

    prop = sol_ds.prop;
    lambda = reshape(prop(:, 1).*prop(:, 2).^2, [1,1,Nz]);
    mu = reshape(prop(:, 1).*prop(:, 3).^2, [1,1,Nz]);

    %%% Evaluate displacement and stress %%%
    % Initialize arrays
    fk_uz  = zeros(Nx,Ny,Nz);  fk_ur  = zeros(Nx,Ny,Nz);
    fk_szz = zeros(Nx,Ny,Nz);  fk_srz = zeros(Nx,Ny,Nz);

    % Wavenumber grid [rad/km]
    [Kx, Ky] = meshgrid(kx, ky);  Kx = Kx';  Ky = Ky';
    Kr = sqrt(Kx.^2 + Ky.^2);     clear Kx Ky;

    for iz = 1:Nz

        % Convert to Cartesian coordinate
        ds1 = sol_ds.ds1(:,:,iz);   ds2 = sol_ds.ds2(:,:,iz);
        ds1_xy = zeros(Nx, Ny, 4);  ds2_xy = zeros(Nx, Ny, 4);
        parfor j = 1:4
            ds1_xy(:,:,j) = interp1(kr, ds1(:,j), Kr, 'linear', 0);
            ds2_xy(:,:,j) = interp1(kr, ds2(:,j), Kr, 'linear', 0);
        end

        % Surface displacement
        fk_uz(:,:,iz) = coeff.c1 .* ds1_xy(:,:,2) + coeff.c2 .* ds2_xy(:,:,2);
        fk_ur(:,:,iz) = coeff.c1 .* ds1_xy(:,:,1) + coeff.c2 .* ds2_xy(:,:,1);
    
        % Surface stress
        if include_stress
            fk_szz(:,:,iz) = coeff.c1 .* ds1_xy(:,:,4) + coeff.c2 .* ds2_xy(:,:,4);
            fk_srz(:,:,iz) = coeff.c1 .* ds1_xy(:,:,3) + coeff.c2 .* ds2_xy(:,:,3);
        end

    end

    %%% Output struct %%%
    % Polar coordinate (NOT FINISHED)
    if strcmp(coord, 'rtz')
        output.uz =  real(ifft(ifft(fk_uz,[],1),[],2)) ./ (dx*dy);
        % output.ur = -real(ifft(ifft(fk_ur,[],1),[],2)) ./ (dx*dy);
        output.szz =  real(ifft(ifft(fk_szz,[],1),[],2)) ./ (dx*dy);
        % output.srz = -real(ifft(ifft(fk_srz,[],1),[],2)) ./ (dx*dy);

    % Cartesian coordinate
    elseif strcmp(coord, 'xyz')

        % Fix k = 0 component
        Kr(1, 1) = Inf;

        % Solve horizontal displacement
        fk_ux = 1j.*fk_ur.*kx./Kr;  fk_uy = 1j.*fk_ur.*ky./Kr;

        output.ux = real(ifft(ifft(fk_ux,[],1),[],2)) ./ (dx*dy);
        output.uy = real(ifft(ifft(fk_uy,[],1),[],2)) ./ (dx*dy);
        output.uz = real(ifft(ifft(fk_uz,[],1),[],2)) ./ (dx*dy);

        if include_stress
            % Solve shear stress
            fk_sxz = 1j.*fk_srz.*kx./Kr;  fk_syz = 1j.*fk_srz.*ky./Kr;
            fk_sxy = mu .* 1j.*(ky.*fk_ux + kx.*fk_uy);
    
            % Solve normal stress
            fk_sxx = (lambda.*fk_szz + 1j.*(4.*mu.*(lambda+mu).*kx.*fk_ux ...
                + 2.*lambda.*mu.*ky.*fk_uy)) ./ (lambda+2.*mu);
            fk_syy = (lambda.*fk_szz + 1j.*(4.*mu.*(lambda+mu).*ky.*fk_uy ...
                + 2.*lambda.*mu.*kx.*fk_ux)) ./ (lambda+2.*mu);

            output.sxx = real(ifft(ifft(fk_sxx,[],1),[],2)) ./ (dx*dy);
            output.syy = real(ifft(ifft(fk_syy,[],1),[],2)) ./ (dx*dy);
            output.szz = real(ifft(ifft(fk_szz,[],1),[],2)) ./ (dx*dy);
            output.sxy = real(ifft(ifft(fk_sxy,[],1),[],2)) ./ (dx*dy);
            output.sxz = real(ifft(ifft(fk_sxz,[],1),[],2)) ./ (dx*dy);
            output.syz = real(ifft(ifft(fk_syz,[],1),[],2)) ./ (dx*dy);
        end

    else
        error('Coordinate system should be rtz or xyz!');
    end
    
    % Spatial axes
    output.x = sol_ds.xh;  output.y = sol_ds.yh;  output.z = sol_ds.zq;

end