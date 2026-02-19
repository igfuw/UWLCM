
### 💧 Microphysics Schemes

Microphysical schemes available in UWLCM are described in detail in the [libcloudphxx documentation](https://github.com/AgnieszkaMakulska/libcloudphxx/tree/docs/docs).

---

### ⚙️ Cases

Cases in UWLCM are pre-configured simulation scenarios. They are implemented as C++ header files in `src/cases/` and include:
- Common functionality through `CasesCommon.hpp`
- Anelastic approximation support through `Anelastic.hpp`
- Various detail implementations in the `detail/` subdirectory

#### **Thermal Cases:**
- **`dry_thermal`** - Dry thermal bubble simulation (rising bubble of dry air)
- **`moist_thermal`** - Moist thermal bubble simulation based on Grabowski & Clark 1999

#### **Atmospheric Boundary Layer Cases:**
- **`dry_pbl`** - Dry planetary boundary layer simulation
- **`dycoms_rf01`** - DYCOMS-II Research Flight 01 (marine stratocumulus case)
- **`dycoms_rf02`** - DYCOMS-II Research Flight 02 (marine stratocumulus case)

#### **Cumulus Congestus Cases:**
- Shared cumulus functionality through `CumulusCongestusCommon.hpp`
- **`cumulus_congestus_icmw20`** - Cumulus congestus simulation (for International Cloud Modeling Workshop 2020)
- **`cumulus_congestus_icmw24`** - Cumulus congestus simulation (for International Cloud Modeling Workshop 2024)

#### **Trade Wind Cumulus Cases:**
- **`rico11`** - RICO (Rain in Cumulus over the Ocean) case
- **`bomex03`** - BOMEX (Barbados Oceanographic and Meteorological Experiment) case
---


### 🧩 Core Components of UWLCM

- **Solvers** (`src/solvers`): numerical algorithms for fluid dynamics and microphysics
- **Forcings** (`src/forcings`): external atmospheric forcing mechanisms
- **Formulae** (`src/formulae`): mathematical formulations for physical processes
- **Detail** (`src/detail`): low-level implementation details and utilities
- **Opts** (`src/opts`): command-line options for simulations

