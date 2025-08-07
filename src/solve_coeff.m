%% Function: Solve coefficients from surface conditions (P-SV system)

function output = solve_coeff(src, sol_ds)

    %%% Read input struct %%%
    Nx = length(sol_ds.xh);  Ny = length(sol_ds.yh);
    dx = sol_ds.dx;  dy = sol_ds.dy;
    kx = sol_ds.kx;  ky = sol_ds.ky;  kr = sol_ds.kr;
    ds1_surf = sol_ds.ds1(:,:,1);  ds2_surf = sol_ds.ds2(:,:,1);


    %%% Convert to Cartesian coordinate %%%
    % Wavenumber grid [rad/km]
    [Kx, Ky] = meshgrid(kx, ky);  Kx = Kx';  Ky = Ky';
    Kr = sqrt(Kx.^2 + Ky.^2);     clear Kx Ky;
    
    % Initialize arrays
    ds1_surf_xy = zeros(Nx, Ny, 4);  ds2_surf_xy = zeros(Nx, Ny, 4);
    
    % Interpolation
    parfor j = 1:4
        ds1_surf_xy(:,:,j) = interp1(kr, ds1_surf(:,j), Kr, 'linear', 0);
        ds2_surf_xy(:,:,j) = interp1(kr, ds2_surf(:,j), Kr, 'linear', 0);
    end


    %%% Solve linear system %%%
    % Normal traction (Positive for tensile)
    fk_pp = -fft(fft(src.pp,[],1),[],2) .* dx*dy;

    % Matrix of the linear system (pagewise)
    A2 = cat(4, ds1_surf_xy(:,:,3:4), ds2_surf_xy(:,:,3:4));

    % Vector of the linear system (pagewise)
    b2 = zeros(1, Nx, Ny, 2);
    b2(1,:,:,:) = cat(3, zeros(Nx, Ny), fk_pp);
    
    % Solve linear system (pagewise)
    c = pagemldivide(permute(A2, [3,4,1,2]), permute(b2, [4,1,2,3]));
    c = squeeze(permute(c, [3,4,1,2]));
    
    % Remove NaN (REMOVE k = 0 component)
    c(1,1,:) = 0;


    %%% Output struct %%%
    output.kx = kx;  output.ky = ky;
    output.c1 = c(:,:,1);  output.c2 = c(:,:,2);
end