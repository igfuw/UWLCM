# University of Warsaw Lagrangian Cloud Model (UWLCM)

🌦️ Welcome to the UWLCM documentation. This is a work in progress.

## Overview
UWLCM is a numerical modeling tool for cloud simulation that combines Large Eddy Simulation (LES) turbulence modeling with Lagrangian cloud microphysics. The model is designed for high-fidelity atmospheric cloud modeling and research applications, enabling researchers to study cloud formation, evolution, and precipitation processes.

## Key Features

- **MPDATA Algorithm**: Advanced finite-difference scheme for accurate advection calculations
- **Different microphysics schemes**:
    - Superdroplet Method for detailed particle-level simulations
    - Bulk microphysics for computational efficiency
- **Turbulence Modeling**:
    - Smagorinsky scheme for explicit SGS modeling
    - Implicit LES for computational efficiency
- **Flexible Dimensionality**: Support for both 2D and 3D simulations
- **GPU Acceleration**: Superdroplet microphysics can utilize GPU computing for enhanced performance



## Dependencies

UWLCM is built upon two key libraries:

- **[libmpdata++](https://github.com/igfuw/libmpdataxx)**: Advanced advection library implementing the MPDATA algorithm
- **[libcloudph++](https://github.com/igfuw/libcloudphxx)**: Cloud microphysics library

## Scientific Publications

The model is described in detail in the following peer-reviewed publications:

- [UWLCM 1.0: A C++ library for atmospheric cloud modeling](https://gmd.copernicus.org/articles/12/2587/2019/)
- [UWLCM 2.0: adaptation of a mixed Eulerian–Lagrangian numerical model for heterogeneous computing clusters](https://gmd.copernicus.org/articles/15/4489/2022/)


## Installation

Detailed build instructions can be found in `Installation_guide.md`.

## Visualization

See UWLCM in action with this visualization of a modeled stratocumulus cloud deck:
[https://youtu.be/EZf260ZKqq0](https://youtu.be/EZf260ZKqq0)


## Contributing

We welcome contributions to UWLCM! Please see the documentation in the `docs/` directory for more information on the codebase and development guidelines.