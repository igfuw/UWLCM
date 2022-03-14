#pragma once

#include "setup.hpp"

namespace detail
{
  // container for vertical profiles
  // CAUTION: new profiles have to be added to both structs and in copy_profiles below
  // TODO: try a different design where it is not necessary ?
  // TODO: store profiles in mem (add an alloc 1d function to libmpdata?) and not in params
  struct profiles_t
  {
    setup::arr_1D_t th_e, p_e, rv_e, rl_e, th_ref, rhod, rhod_vctr, w_LS, hgt_fctr, th_LS, rv_LS, mix_len, relax_th_rv_coeff;
    std::array<setup::arr_1D_t, 2> geostr;
  
    profiles_t(int nz) :
    // rhod_vctr is bigger, because it is density at the position of vector fields (Arakawa-C)
      th_e(nz), p_e(nz), rv_e(nz), rl_e(nz), th_ref(nz), rhod(nz), rhod_vctr(nz+1), w_LS(nz), hgt_fctr(nz), th_LS(nz), rv_LS(nz), mix_len(nz), relax_th_rv_coeff(nz)
    {
      geostr[0].resize(nz);
      geostr[1].resize(nz);
  
      // set to zero just to have predicatble output in cases that dont need these profiles
      geostr[0]     = 0.;
      geostr[1]     = 0.;
      hgt_fctr      = 0.;
      rl_e          = 0.;
      relax_th_rv_coeff = 0.;
    }
  };
  
  struct profile_ptrs_t
  {
    setup::arr_1D_t *th_e, *p_e, *rv_e, *rl_e, *th_ref, *rhod, *rhod_vctr, *w_LS, *hgt_fctr, *geostr[2], *th_LS, *rv_LS, *mix_len, *relax_th_rv_coeff;
  };
  
  // copy external profiles into rt_parameters
  // TODO: more elegant way
  template<class params_t>
  inline void copy_profiles(profiles_t &profs, params_t &p)
  {
    std::vector<std::pair<std::reference_wrapper<setup::arr_1D_t*>, std::reference_wrapper<setup::arr_1D_t>>> tobecopied = {
      {p.hgt_fctr          , profs.hgt_fctr          },
      {p.th_e              , profs.th_e              },
      {p.p_e               , profs.p_e               },
      {p.rv_e              , profs.rv_e              },
      {p.rl_e              , profs.rl_e              },
      {p.th_ref            , profs.th_ref            },
      {p.rhod              , profs.rhod              },
      {p.rhod_vctr         , profs.rhod_vctr         },
      {p.w_LS              , profs.w_LS              },
      {p.th_LS             , profs.th_LS             },
      {p.rv_LS             , profs.rv_LS             },
      {p.geostr[0]         , profs.geostr[0]         },
      {p.geostr[1]         , profs.geostr[1]         },
      {p.mix_len           , profs.mix_len           },
      {p.relax_th_rv_coeff , profs.relax_th_rv_coeff }
    };
  
    for (auto dst_src : tobecopied)
    {
      dst_src.first.get() = new setup::arr_1D_t(dst_src.second.get().dataFirst(), dst_src.second.get().shape(), blitz::neverDeleteData);
    }
  }
};
