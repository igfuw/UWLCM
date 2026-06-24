### How to create a new simulation case?


To define a new physical scenario, a new file must be added to `src/cases/`. 
Initial and boundary conditions, surface fluxes and domain size should be defined there.
The choice of the new case must be enabled in `src/run_hplr.cpp`.
When running the simulation, one of the cases is selected as a command line parameter.


### Handling data types and units

- Real numbers are represented by the `real_t` variable. It enables the choice between different precision types (float or double).
- Physical quantities are represented by `quantity<si::variable_type, real_t>` (defined in the Boost C++ library). 
`variable_type` can be one of  `dimensionless`, `length`, `area`, `mass`, `time`, `temperature`, and other SI quantities. Quantities such as 1/m can be created with the variable type:
`divide_typeof_helper<si::dimensionless, si::length>::type`.

