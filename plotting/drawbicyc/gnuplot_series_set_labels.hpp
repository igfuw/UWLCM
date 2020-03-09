#include "common.hpp"

void gnuplot_series_set_labels(Gnuplot &gp, std::string plt)
{
  gp << "set yrange[*:*]\n";
  gp << "set xrange[*:*]\n";

  if (plt == "clfrac")
  {
  //  res_pos *= 60.;
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set title 'cloud fraction'\n";
  }
  else if (plt == "ract_com")
  {
    gp << "set ylabel 'q_c center of mass [km]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'center of mass of cloud water'\n";
  }
  else if (plt == "th_com")
  {
    gp << "set ylabel 'th prtrb center of mass [km]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'center of mass height'\n";
  }
  else if (plt == "ract_avg")
  {
    gp << "set ylabel '<q_c>  [g/kg]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'avg cloud water mixing ratio'\n";
  }
  else if (plt == "ract_std_dev")
  {
    gp << "set ylabel 'sigma(q_c) / <q_c>'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'relative std dev of q_c'\n";
  }
  else if (plt == "cloud_avg_act_conc")
  {
    gp << "set ylabel '<N_c>  [1/cm^3]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'avg concentration of activated droplets'\n";
  }
  else if (plt == "cloud_std_dev_act_conc")
  {
    gp << "set ylabel 'sigma(N_c) / <N_c>'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'relative std dev N_c'\n";
  }
  else if (plt == "cloud_avg_supersat")
  {
    gp << "set ylabel '<S> [%]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'avg supersaturation'\n";
  }
  else if (plt == "cloud_std_dev_supersat")
  {
    gp << "set ylabel 'sigma(S) / <S>'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'relative std dev S'\n";
  }
  else if (plt == "cloud_avg_act_rad")
  {
    gp << "set ylabel '<r_{mean}> [um]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'average mean radius of activated droplets'\n";
  }
  else if (plt == "cloud_avg_std_dev_act_rad")
  {
    gp << "set ylabel '<sigma(r)> [um]'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'average std dev of radius of activated droplets'\n";
  }
  else if (plt == "sd_conc_avg")
  {
    gp << "set ylabel '<N_{SD}>'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'average SD number'\n";
  }
/*
    else if (plt == "sd_conc_std_dev")
    {
      res_pos *= 60.;
      gp << "set ylabel 'sigma(N_{SD}) / <N_{SD}>'\n";
      gp << "set xlabel 'time [min]'\n";
      gp << "set title 'relative std dev of N_{SD}'\n";
    }
*/
  else if (plt == "sd_conc_act_avg")
  {
    gp << "set ylabel '<N_{SD}^{act}>'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set title 'average activated SD number'\n";
  }
  else if (plt == "tot_water")
  {
    gp << "set title 'mean water mass density'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'mean (hydrometeors+vapor) water mass density [g/m^3]'\n";
  }
  else if (plt == "com_vel")
  {
    gp << "set title 'vertical velocity of the COM\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'w [m/s]'\n";
  }
  else if (plt == "com_supersat")
  {
    gp << "set title 'supersaturation at the COM\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'S [%]'\n";
  }
  else if (plt == "com_mom0")
  {
    gp << "set title 'cloud drops concentration at the center of mass\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'N_c [cm^{-3}]'\n";
  }
  else if (plt == "com_mom1")
  {
    gp << "set title 'mean wet radius at the center of mass\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel '<r> [micron]'\n";
  }
  else if (plt == "com_mom2")
  {
    gp << "set title 'std dev of r at the center of mass\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'std dev [micron]'\n";
  }
  else if (plt == "com_sd_conc")
  {
    gp << "set title 'number of SDs at the center of mass\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'SD no'\n";
  }
  else if (plt == "nc")
    gp << "set title 'average cloud drop conc [1/cm^3]'\n";
  else if (plt == "cl_nc")
  {
    gp << "set title 'average cloud drop conc [1/cm^3] in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_nr")
  {
    gp << "set title 'average rain drop conc [1/cm^3] in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "wvarmax")
  {
    gp << "set title 'max variance of w [m^2 / s^2]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "surf_precip")
  {
    gp << "set title 'surface precipitation [mm/d]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "acc_precip")
  {
    gp << "set title 'accumulated surface precipitation [mm]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "mass_dry")
    gp << "set title 'total dry mass [g]'\n";
  else if (plt == "RH_max")
  {
    gp << "set title 'max RH in the domain'\n";
    gp << "set xlabel 'time [min]'\n";
    gp << "set ylabel 'RH'\n";
  }
  else if (plt == "lwp")
  {
    gp << "set title 'liquid water path [g / m^2]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "rwp")
  {
    gp << "set title 'rain water path [g / m^2]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cloud_base")
  {
    gp << "set title 'cloud base [m]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_gccn_conc")
  {
    gp << "set title 'average gccn conc [1/cm^3] in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_gccn_to_non_gccn_conc_ratio")
  {
    gp << "set title 'average ratio of gccn to non-gccn conc in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_non_gccn_conc")
  {
    gp << "set title 'average non gccn conc [1/cm^3] in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_gccn_meanr")
  {
    gp << "set title 'average wet radius [um] of GCCNs in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_meanr")
  {
    gp << "set title 'average wet radius [um] of cloud droplets in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "cl_avg_cloud_rad")
  {
    gp << "set title 'average wet radius [um] of cloud droplets in cloudy cells'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "gccn_conc")
  {
    gp << "set title 'average gccn conc [1/cm^3]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "non_gccn_conc")
  {
    gp << "set title 'average non gccn conc [1/cm^3]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "surf_flux_sensible")
  {
    gp << "set title 'sensible surf flux [W/m^2]'\n";
  }
  else if (plt == "surf_flux_latent")
  {
    gp << "set title 'latent surf flux [W/m^2]'\n";
  }
  else if (plt == "er")
  {
    gp << "set title 'entrainment rate [cm / s]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "tot_tke")
  {
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
    gp << "set title 'turbulent kinetic energy (resolved + sgs) [m^3 / s^2]'\n";
  }
  else if (plt == "clb_bigrain_mean_inclt")
  {
    gp << "set title 'big rain (r>40um) at clbase: <time since activation> [s]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "clb_bigrain_mean_rd")
  {
    gp << "set title 'big rain (r>40um) at clbase: <dry radius> [um]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "clb_bigrain_mean_kappa")
  {
    gp << "set title 'big rain (r>40um) at clbase: <kappa>'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "clb_bigrain_mean_conc")
  {
    gp << "set title 'big rain (r>40um) at clbase: <conc.> [1/cc]'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
  else if (plt == "clb_bigrain_mean_gccn_fraction")
  {
    gp << "set title 'big rain (r>40um) at clbase: <fraction with kappa>0.61>'\n";
    gp << "set xlabel ''\n";
    gp << "set ylabel ''\n";
  }
}
