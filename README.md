# University of Warsaw Lagrangian Cloud Model

UWLCM is a tool for numerical modeling of clouds using LES model of turbulence and Lagrangian cloud microphysics.
Paper describing the model is available at: https://www.geosci-model-dev-discuss.net/gmd-2018-281/

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
 
Visualization of a modeled stratocumulus cloud deck:
https://youtu.be/EZf260ZKqq0
