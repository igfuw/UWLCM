# University of Warsaw Lagrangian Cloud Model

UWLCM is a tool for numerical modeling of clouds using LES model of turbulence and Lagrangian cloud microphysics.

Papers describing the model: 
https://gmd.copernicus.org/articles/12/2587/2019/
https://gmd.copernicus.org/articles/15/4489/2022/


Key features:
 - MPDATA algorithm for advection
 - Superdroplet or bulk microphysics
 - SGS turbulence modeled with Smagorinsky scheme or implicit LES
 - models for SGS motion, condensation and coalescence of superdroplets
 - 2D or 3D simulations
 - superdroplet microphysics can be calculated using GPUs

UWLCM is implemented using the libmpdata++ advection library and libcloudph++ microphysics library.
https://github.com/igfuw/libmpdataxx
https://github.com/igfuw/libcloudphxx

Build instructions can be found in build_instructions.rtf.
 
Visualization of a modeled stratocumulus cloud deck:
https://youtu.be/EZf260ZKqq0
