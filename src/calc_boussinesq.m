%% Function: Boussinesq solution of displacement & stress

function output = calc_boussinesq(mesh, elast_prop, coord)

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

    % Polar distance [km]
    r_2d = sqrt((xh-x0).^2 + (yh-y0).^2);
    r_3d = sqrt(r_2d.^2 + zq.^2);


    %%% Boussinesq solution %%%
    % Coordinate system: +z-direction upward
    % Stress convention: Positive for tensile stress

    % Displacement [μm]
    ur =  1/(4*pi*mu) .* r_2d .* (zq./r_3d.^3 - (1-2*nu)./r_3d./(zq+r_3d));
    uz = -1/(4*pi*mu) .* (2*(1-nu)./r_3d + zq.^2./r_3d.^3);

    % Stress [Pa]
    srr = 1/(2*pi) .* ((1-2*nu)./r_3d./(zq+r_3d) - 3.*r_2d.^2.*zq./r_3d.^5);
    stt = (1-2*nu)/(2*pi) .* (zq./r_3d.^3 - 1./r_3d./(zq+r_3d));
    szz = -3/(2*pi) .* zq.^3./r_3d.^5;
    srz =  3/(2*pi) .* r_2d.*zq.^2./r_3d.^5;
    

    %%% Output struct %%%
    % Polar coordinate
    if strcmp(coord, 'rtz')
        output.ur = ur;     output.uz = uz;
        output.srr = srr;   output.stt = stt;
        output.szz = szz;   output.srz = srz;

    % Cartesian coordinate
    elseif strcmp(coord, 'xyz')

        % Convert displacement components
        cs = (xh-x0) ./ r_2d;   ss = (yh-y0) ./ r_2d;
        ux = ur .* cs;          uy = ur .* ss;

        % Convert stress components
        sxx = srr.*cs.^2 + stt.*ss.^2;  syy = stt.*cs.^2 + srr.*ss.^2;
        sxz = srz.*cs;  syz = srz.*ss;  sxy = (srr-stt).*cs.*ss;

        output.uz = uz;     output.ux = ux;     output.uy = uy;
        output.sxx = sxx;   output.syy = syy;   output.szz = szz;
        output.sxz = sxz;   output.syz = syz;   output.sxy = sxy;

    else
        error('Coordinate system should be rtz or xyz!');
    end

    % Spatial axes
    output.x = xh;  output.y = yh;  output.z = zq;

end