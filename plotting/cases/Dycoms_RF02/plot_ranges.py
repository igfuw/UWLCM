xscaledict = {
  "thl" : "linear",
  "00rtot" : "linear",
  "rliq" : "linear",
  "prflux" : "linear",
  "cl_nc" : "linear",
  "clfrac" : "linear",
  "wvar" : "linear",
  "w3rd" : "linear",
  "sat_RH" : "linear",
  "rad_flx" : "linear",
  "lwp" : "linear",
  "er" : "linear",
  "wvarmax" : "linear",
  "surf_precip" : "linear",
  "acc_precip" : "linear",
  "cloud_base" : "linear",
  "gccn_rw_cl" : "linear",
  "non_gccn_rw_cl" : "linear",
  "base_prflux_vs_clhght" : "linear",
  "cl_gccn_conc" : "linear"
}

yscaledict = {
  "thl" : "linear",
  "00rtot" : "linear",
  "rliq" : "linear",
  "prflux" : "linear",
  "cl_nc" : "linear",
  "clfrac" : "linear",
  "wvar" : "linear",
  "w3rd" : "linear",
  "sat_RH" : "linear",
  "rad_flx" : "linear",
  "lwp" : "linear",
  "er" : "linear",
  "wvarmax" : "linear",
  "surf_precip" : "linear",
  "acc_precip" : "linear",
  "cloud_base" : "linear",
  "gccn_rw_cl" : "linear",
  "non_gccn_rw_cl" : "linear",
  "base_prflux_vs_clhght" : "log",
  "cl_gccn_conc" : "log"
}

xlimdict_profs = {
  "thl" : None,#(288.2,289.2),
  "00rtot" : None,#(9,10.4),
  "rliq" : (0,0.7),
  "prflux" : None,#(0,70),
  "cl_nc" : None,#(0,120),
  "clfrac" : None,
  "wvar" : (-0.1,0.8),
  "w3rd" : (-0.15,.15),
  "sat_RH" : (-5,1),
  "rad_flx" : None,
  "gccn_rw_cl" : None,#(0,40),
  "non_gccn_rw_cl" :None,# (0,7),
  "base_prflux_vs_clhght" : None,#(1,10000)
}

ylimdict_profs = {
  "thl" : (0,1.2),
  "00rtot" : (0,1.2),
  "rliq" : (0,1.2),
  "prflux" : (0,1.2),
  "cl_nc" : (0,1.2),
  "clfrac" : (0,1.2),
  "wvar" : (0,1.2),
  "w3rd" : (0,1.2),
  "sat_RH" : (0,1.2),
  "rad_flx" : (0,1.2),
  "gccn_rw_cl" : (0,1.2),
  "non_gccn_rw_cl" : (0,1.2),
  "base_prflux_vs_clhght" : (0,2500)
}

xlimdict_series = {
  "clfrac" : (1,5),
  "cl_nc" : (1,5),
  "lwp" : (1,5),
  "er" : (1,5),
  "wvarmax" : (1,5),
  "surf_precip" : (1,5),
  "acc_precip" : (1,5),
  "cloud_base" : (1,5),
  "cl_gccn_conc" : (1,5)
}

ylimdict_series = {
  "clfrac" : None,
  "cl_nc" : None,
  "lwp" : None,
  "er" : None,
  "wvarmax" : None,
  "surf_precip" : None,#(-0.01,0.5),
  "acc_precip" : None,#(0,0.07),
  "cloud_base" : None,
  "cl_gccn_conc" : None,#(1e-10,1e-0)
}


