# Quasi-Static Deformation Modeling

Modeling the quasi-static deformation of a layered elastic medium under surface loading using the **Propagator Matrix Method**.

**Author:** Qing Ji, Stanford University  (qingji@stanford.edu)

---

## Overview
This repository contains MATLAB scripts and functions for modeling the deformation of a layered elastic medium under static or quasi-static loading at the surface. The solution is obtained via the **Propagator Matrix Method**.

---

## Requirements
- MATLAB R202x or later
- [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html) *(optional)*  
  You can replace all `parfor` loops with `for` loops in:
  - `solve_ds.m`
  - `calc_layer.m`

---

## Repository Structure

### 1. Scripts
- **`main_example.m`** — Example workflow that illustrates the usage of functions
- **`test_Boussinesq.m`, `test_Sorrells.m`** — Benchmark for Boussinesq and Sorrells problems (point load and pressure wave applied at the surface of elastic halfspace)

### 2. MATLAB Functions (`./src`)
- **`solve_ds.m`**, **`calc_layer.m`** — Core propagator matrix method
- **`create_psfc.m`** — Example surface pressure loading (Gaussian, pressure wave, turbulent pressure, Delta load)

### 3. Data (`./data`)
- **`rwb_cb.mat`** — Red-white-blue colorbar  
- **`vel_model.csv`** — Example layered medium over halfspace  
- **`cm1out.nc`** — Turbulent pressure field from large eddy simulation  
  *(Download: [Google Drive link](https://drive.google.com/file/d/19DOqOyQnbwYKHG0U1_NNDzsCQnU6Bwfq/view?usp=sharing))*

---

## Notes

The propagator matrix method uses a spatial horizontal Fourier transform, making it ideal for **periodic loading**.  
For concentrated loads, set the spatial domain much larger than the loaded area.

The solution neglects the zero-wavenumber component (`k = 0`). Details of the propagator matrix method can be found in [Method.pdf](./Method.pdf).
