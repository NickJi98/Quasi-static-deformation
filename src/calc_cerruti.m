%% Function: Cerruti solution of displacement & stress

function output = calc_cerruti(mesh, elast_prop, coord)

    %%% Elastic properties %%%
    rho = elast_prop(1);  vp = elast_prop(2);  vs = elast_prop(3);

    % Shear modulus [GPa] & Poisson's ratio
    mu = rho * vs^2;  nu = ((vp/vs)^2-2) / ((vp/vs)^2-1) / 2;


    %%% Mesh grid %%%
    % Grid points on horizontal axis [km]
    xh = (0:mesh.Nx-1)' .* mesh.dx;  x0 = mesh.Nx * mesh.dx / 2;
    try
        yh = (0:mesh.Ny-1) .* mesh.dy;  y0 = mesh.Ny * mesh.dy / 2;
    catch
        yh = xh';  y0 = x0;
    end

    % Depth query points [km]
    zq = (0:mesh.Nz-1)' .* mesh.dz;  zq = reshape(zq, [1,1,mesh.Nz]);

    % Relative position & distance [km]
    xr = xh - x0;  yr = yh - y0;
    r_3d = sqrt(xr.^2 + yr.^2 + zq.^2);


    %%% Cerruti solution %%%
    % Vertical displacement: Positive for upward motion
    % Stress convention: Positive for tensile stress
    % Horizontal load Q = 1e6 N along +x direction

    % Stresses are derived from the displacement field with Hooke's law
    % (formulas verified symbolically: equilibrium & traction-free surface)

    % Displacement [μm]
    ux =  1/(4*pi*mu) .* (1./r_3d + xr.^2./r_3d.^3 ...
        + (1-2*nu).*(1./(r_3d+zq) - xr.^2./r_3d./(r_3d+zq).^2));
    uy =  1/(4*pi*mu) .* xr.*yr ...
        .* (1./r_3d.^3 - (1-2*nu)./r_3d./(r_3d+zq).^2);
    uz = -1/(4*pi*mu) .* xr ...
        .* (zq./r_3d.^3 + (1-2*nu)./r_3d./(r_3d+zq));

    % Stress [Pa]
    sxx = xr./(2*pi.*r_3d.^3) .* (-3.*xr.^2./r_3d.^2 + (1-2*nu)./(r_3d+zq).^2 ...
        .* (r_3d.^2 - yr.^2 - 2.*r_3d.*yr.^2./(r_3d+zq)));
    syy = xr./(2*pi.*r_3d.^3) .* (-3.*yr.^2./r_3d.^2 + (1-2*nu)./(r_3d+zq).^2 ...
        .* (3.*r_3d.^2 - xr.^2 - 2.*r_3d.*xr.^2./(r_3d+zq)));
    sxy = yr./(2*pi.*r_3d.^3) .* (-3.*xr.^2./r_3d.^2 + (1-2*nu)./(r_3d+zq).^2 ...
        .* (-r_3d.^2 + xr.^2 + 2.*r_3d.*xr.^2./(r_3d+zq)));
    szz = -3.*xr.*zq.^2 ./ (2*pi.*r_3d.^5);
    sxz =  3.*xr.^2.*zq ./ (2*pi.*r_3d.^5);
    syz =  3.*xr.*yr.*zq ./ (2*pi.*r_3d.^5);


    %%% Output struct %%%
    % Cartesian coordinate
    if strcmp(coord, 'xyz')
        output.uz = uz;     output.ux = ux;     output.uy = uy;
        output.sxx = sxx;   output.syy = syy;   output.szz = szz;
        output.sxz = sxz;   output.syz = syz;   output.sxy = sxy;

    else
        error('Coordinate system should be xyz for Cerruti solution!');
    end

    % Spatial axes
    output.x = xh;  output.y = yh;  output.z = zq;

end
