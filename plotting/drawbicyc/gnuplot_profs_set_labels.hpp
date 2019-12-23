#include "common.hpp"

void gnuplot_profs_set_labels(Gnuplot &gp, std::string plt, const bool normalize)
{
  gp << "set xrange [*:*]\n";
  if(normalize)
    gp << "set yrange [0.:1.2]\n";
  else
    gp << "set yrange [*:*]\n";

  if (plt == "rliq")
  {
    gp << "set title 'liquid water [g/kg]'\n";
  }
  if (plt == "gccn_rw")
  {
    gp << "set title 'gccn-based droplets mean wet radius'\n";
  }
  if (plt == "non_gccn_rw")
  {
    gp << "set title 'non-gccn-based droplets mean wet radius'\n";
  }
  if (plt == "non_gccn_rw_down")
  {
    gp << "set title 'non-gccn-based droplets mean wet radius (downdraughts only)'\n";
  }
  if (plt == "gccn_rw_down")
  {
    gp << "set title 'gccn-based droplets mean wet radius (downdraughts only)'\n";
  }
  if (plt == "non_gccn_rw_up")
  {
    gp << "set title 'non-gccn-based droplets mean wet radius (updraughts only)'\n";
  }
  if (plt == "gccn_rw_up")
  {
    gp << "set title 'gccn-based droplets mean wet radius (updraughts only)'\n";
  }
  if (plt == "ugccn_rw_down")
  {
    gp << "set title 'ultra-gccn-based droplets mean wet radius (downdraughts only)'\n";
  }
  if (plt == "act_conc_up")
  {
    gp << "set title 'activated droplets concentation [1/cm^3] (updrafts)'\n";
  }
  if (plt == "nc_up")
  {
    gp << "set title 'clloud droplets concentation [1/cm^3] (updrafts)'\n";
  }
  if (plt == "cl_nc_up")
  {
    gp << "set title 'cloud droplets concentation [1/cm^3] (cloudy updrafts)'\n";
  }
  if (plt == "nc_down")
  {
    gp << "set title 'clloud droplets concentation [1/cm^3] (updrafts)'\n";
  }
  if (plt == "act_rd_up")
  {
    gp << "set title 'activated droplets mean dry radius (updrafts)'\n";
  }
  if (plt == "actRH_rd")
  {
    gp << "set title 'activated (RH>Sc) droplets mean dry radius'\n";
  }
  if (plt == "actrw_rd")
  {
    gp << "set title 'activated (rw>rc) droplets mean dry radius'\n";
  }
  else if (plt == "rv")
  {
    gp << "set title 'rv [g/kg] averaged over 2h-6h, w/o rw<0.5um'\n";
    if(normalize)
    {
      gp << "set yrange [0.:0.6]\n";
      gp << "set xrange [9.2:9.5]\n";
    }
  }
  else if (plt == "u")
  {
    gp << "set title 'u [m/s]'\n";
  }
  else if (plt == "v")
  {
    gp << "set title 'v [m/s]'\n";
  }
  else if (plt == "w")
  {
    gp << "set title 'w [m/s]'\n";
  }
  else if (plt == "sd_conc")
  {
    gp << "set title '# of SD'\n";
  }
  else if (plt == "vel_div")
  {
    gp << "set title 'vel_div [1/s]'\n";
  }
  else if (plt == "sat_RH")
  {
    gp << "set title 'supersaturation RH-based'\n";
  }
  else if (plt == "sat_RH_up")
  {
    gp << "set title 'supersaturation RH-based in updrafts only'\n";
    gp << "set yrange [0.45:1.]\n";
      gp << "set xrange [0.000:*]\n";
  }
  else if (plt == "00rtot")
  {
    gp << "set title 'total water [g/kg]'\n";
  }
  else if (plt == "N_c")
  {
    gp << "set title 'cloud droplets ( 0.5um < r < 25um) concentration [1/cm^3]'\n";
  }
  else if (plt == "cl_nc")
  {
    gp << "set title 'cloud droplets concentration in cloudy cells [1/cm^3]'\n";
  }
  else if (plt == "thl")
  {
    gp << "set title 'liquid potential temp [K]'\n";
  }
  else if (plt == "clfrac")
  {
    gp << "set title 'cloud fraction'\n";
  }
  else if (plt == "base_prflux_vs_clhght")
  {
    gp << "set title 'cl base pr flux vs cl hgt'\n";
  }
  else if (plt == "prflux")
  {
    gp << "set title 'precipitation flux [W/m^2]'\n";
  }
  else if (plt == "rad_flx")
  {
    gp << "set title 'radiative flux [W/m2]'\n";
  }
  else if (plt == "wvar")
  {
    gp << "set title 'variance of w [m^2 / s^2]'\n";
  }
  else if (plt == "w3rd")
  {
    gp << "set title '3rd mom of w [m^3 / s^3]'\n";
  }
  else if (plt == "sgs_tke")
  {
    gp << "set title 'sgs tke [TODO]'\n";
  }
  else if (plt == "k_m")
  {
    gp << "set title 'k_m [TODO]'\n";
  }
  else if (plt == "sgs_tht_flux")
  {
    gp << "set title 'sgs tht flux [TODO]'\n";
  }
  else if (plt == "sgs_rv_flux")
  {
    gp << "set title 'sgs rv flux [TODO]'\n";
  }
  else if (plt == "sgs_rc_flux")
  {
    gp << "set title 'sgs rc flux [TODO]'\n";
  }
  else if (plt == "sgs_u_flux")
  {
    gp << "set title 'sgs u flux [TODO]'\n";
  }
}
