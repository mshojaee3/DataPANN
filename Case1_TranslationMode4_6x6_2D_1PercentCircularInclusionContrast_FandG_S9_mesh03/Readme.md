# Hyperelastic Translation Ensemble Dataset

## Problem Statement

This dataset is generated from **2D hyperelastic simulations** of a heterogeneous microstructured material using Abaqus and Python drivers.

The objective is to learn the constitutive mapping:

$$
\left(\mathbf{F}, \nabla\mathbf{F}\right)
\equiv (\mathbf{F}, \mathbf{G})
\longmapsto
\left(\mathbf{P}, \mathbf{Q}\right)
$$

where macroscopic deformation quantities \((F, G)\) are mapped to macroscopic stress responses \((P, Q)\).

The dataset captures:

- Microstructure heterogeneity
- Inclusion translation within a periodic unit cell
- Higher-order deformation and stress effects
- Data-driven constitutive modeling behavior

Each load case is a point in a reduced **9D boundary-condition coefficient space** \((H, G)\), sampled on the unit sphere \(S^8\) using energy-minimizing distributions.

## Geometry

### Specimen

- Total domain size:
  - `Lx_tot = 0.3`
  - `Ly_tot = 0.3`
- Number of cells:
  - `Nx = 6`
  - `Ny = 6`

### Unit Cell

Each unit cell contains a circular inclusion. Its geometric parameters are:

$$
R = 0.30\,L_{\mathrm{cell}},
\qquad
h = 0.03\,L_{\mathrm{cell}},
\qquad
\ell_{\mathrm{ligament}} = 0.10\,L_{\mathrm{cell}}
$$

### Material Contrast

- Heterogeneous material: `kappaMat = 0.01`
- Homogeneous reference material: `kappaMat = 1.0`

## Translation Ensemble

The dataset supports translations of the inclusion within each unit cell to represent microstructural variability.

The inclusion-center translation is denoted by:

$$
\boldsymbol{\zeta} = (\zeta_1,\zeta_2)
$$

Its admissible magnitude is bounded by:

$$
\zeta_{\max} = \frac{L}{2} - R - m,
$$

where \(m\) is a prescribed boundary margin.

Supported translation modes:

1. Diagonal
2. Full grid
3. Cross
4. Center only

The current setup uses **center-only mode**, meaning inclusions remain centered and no translation variation is included.

## Load-Case Space

The deformation is parameterized by:

- `H`: linear deformation terms
- `G`: quadratic deformation terms

The original 12D parameter space is reduced to **9D** because:

1. \(G_{ijk}\) is symmetric in \(j,k\).
2. The skew-symmetric part of \(H\) represents rigid rotation and does not contribute to energy.

The retained coefficient vector is:

$$
\mathbf{c} =
\left[
H_{11},\; H_{22},\; H_{12},\;
G_{1,11},\; G_{1,22},\; G_{1,12},\;
G_{2,11},\; G_{2,22},\; G_{2,12}
\right]^{\mathsf{T}}
\in \mathbb{R}^{9}.
$$

Total dimension: **9**.

## Simulation Workflow

Each load case is simulated twice.

### Heterogeneous Simulation

Uses the actual microstructure with `kappaMat != 1` and produces stress, energy, reaction forces, and deformation measures.

### Homogeneous Simulation

Uses the same geometry with uniform material properties. It provides a smooth deformation field for pairing and input consistency.

### Pairing

- Frames are matched exactly in time.
- Only valid heterogeneous-homogeneous frame pairs are retained.
- Final CSV values are computed from the heterogeneous simulation.

## Data Structure

Each load case produces:

```text
TrainingData/<LoadCaseName>_FINAL.csv
```

Each row corresponds to one simulation time frame.

### Deformation Features

Input features:

- `F11`, `F12`, `F21`, `F22`: mean deformation gradient.
- `G111` through `G222`: eight components describing the spatial gradient of \(F\).

### Stress Targets

- `P11`, `P12`, `P21`, `P22`: first Piola-Kirchhoff stress.
- `Q111` through `Q222`: eight components of the first spatial moment of stress.

### Energy

- `ALLSE_mean`: total strain energy
- `W_mean`: energy density

### Diagnostics

- `Jmin`: minimum determinant of \(F\)
- `FfitResid`: affine approximation residual
- `Area`: domain area

### Reaction Forces

Eight reaction-force components:

- Left boundary: X/Y
- Right boundary: X/Y
- Bottom boundary: X/Y
- Top boundary: X/Y

### Metadata

Constant for each load case:

- `NX`, `NY`
- `Lx_tot`, `Ly_tot`
- `Rfrac`, `meshFrac`
- `kappaMat`
- `MaterialID`
- `LoadCase`
- `NFRAMES`

## Physical Quantities

### Deformation Gradient

$$
\mathbf{F} = \mathbf{I} + \nabla \mathbf{u}
$$

### First Piola-Kirchhoff Stress

$$
\mathbf{P} = J\,\boldsymbol{\sigma}\,\mathbf{F}^{-\mathsf{T}},
\qquad
J = \det(\mathbf{F})
$$

### First Spatial Moment of Stress

$$
Q_{ikl} = \left\langle P_{ik}\,X_c^l \right\rangle
$$

### Energy Density

$$
W = \frac{\mathrm{ALLSE}}{V_0}
$$

## Usage

Typical applications include:

- Training neural-network constitutive models: input \((F, G)\), output \((P, Q)\)
- Nonlinear homogenization studies
- Higher-order continuum modeling
- Microstructure sensitivity analysis
- Translation-invariance studies

## Notes

- Load cases are hierarchical: increasing the dataset size preserves earlier cases.
- Frame consistency is strictly enforced; mismatched frames are discarded.
- Partial simulations may be included if they satisfy the minimum-frame threshold.
- All final deformation, stress, energy, and reaction-force values are taken from heterogeneous simulations.
- Homogeneous simulations are used only for pairing and consistency checks.

## Summary

This dataset provides a physics-consistent mapping:

$$
\left(\mathbf{F}, \mathbf{G}\right)
\longmapsto
\left(\mathbf{P}, \mathbf{Q}\right)
$$

for heterogeneous hyperelastic microstructures, including higher-order effects and an extensible inclusion-translation ensemble.
