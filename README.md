# Quasi-static-deformation
Modeling the quasi-static deformation of a layered elastic medium under surface loading using Propagator Matrix Method

Author: Qing Ji, Stanford University (qingji@stanford.edu)

Still updating ...


MATLAB packages:

I use MATLAB Parallel Computing Toolbox to convert loops into the parallel version (i.e., 'parfor' instead of 'for' loop). You may install this toolbox, or just simply replace all 'parfor' loop with 'for' loop in functions solve_ds and solve_coeff



1. Scripts

1.1 test_halfspace_xyz.m: Benchmark with Boussinesq solution (i.e., point load over an elastic halfspace)

1.2 main_example.m: Example workflow and inputs



2. MATLAB functions (./my_func)

2.1 calc_boussinesq: Boussinesq solution (see ./docs/Boussinesq.pdf for expressions; note the difference in configuration)

2.2 solve_ds, solve_coeff, calc_layer: Propagator matrix method (see ./docs/Propagator.pdf, section 5.5 & 5.6)

2.3 create_psfc: Example surface pressure loading (Gaussian, single Fourier mode, turbulent pressure, Delta load from white spectrum) 

2.4 plot_2d, plot_compare_2d, plot_result_2d: Plotting functions



3. Notes

- 3.1 Since propagator matrix method is based on spatial horizontal Fourier transform, periodic loading is best represented. 

However, for a concentrated loading, you may set the spatial domain much larger than the dominantly loaded area.

- 3.2 The zero-wavenumber (k = 0, DC) component is neglected in the propagator matrix method. This is because a uniform loading over an infinite halfspace just leads to an offset of the top surface, and this part is neglected. 

In short, the numerical result will have zero means. For benchmark with the Boussinesq solution, I add the mean of the analytical solution in the domain to the numerical result.

- 3.3 Numerical stability. The propagator matrix method starts from the bottom halfspace, and use matrix exponential (i.e., exp(Ah) where A is a matrix and h is the layer thickness) to propagate upward until reaching the surface. 

Therefore, for larger wavenumber k, the 'stability depth' becomes shorter, so you don't want the layered medium to be too thick.

In short, for small k (large wavelength) you can make the layered medium more coarse and deeper.



4. Data (./my_data)

4.1 rwb_cb.mat: Colorbar red-white-blue

4.2 vel_model.csv: Example layered medium over a bottom halfspace

4.3 cm1out.nc: Turbulent pressure field from large eddy simulation (download here: https://drive.google.com/file/d/19DOqOyQnbwYKHG0U1_NNDzsCQnU6Bwfq/view?usp=sharing)