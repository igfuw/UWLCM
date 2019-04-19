#pragma once

const std::vector<std::string> series_dycoms({
"wvarmax", "clfrac", "lwp", "er",
 "surf_precip", 
//"mass_dry", 
 "acc_precip",
 "cl_nc",
 "cloud_base"
// "sd_conc_avg", "sd_conc_std_dev",
// "tot_water"
});

const std::vector<std::string> series_moist_thermal({
"cloud_avg_act_conc", //"cloud_std_dev_act_conc",
"ract_avg", //"ract_std_dev",
"sd_conc_avg", //"sd_conc_std_dev",
"sd_conc_act_avg", //"sd_conc_std_dev",
"cloud_avg_supersat",//"cloud_std_dev_supersat",
"cloud_avg_act_rad",
"cloud_avg_std_dev_act_rad",
"RH_max",
"ract_com",
"clfrac",
 "tot_water",
"com_mom0","com_mom1","com_mom2", // higher moments need lower ones enabled!!
"com_vel",
"com_supersat",
"com_sd_conc"

// "mass_dry",
//"th_com"
});

const std::vector<std::string> series_sgs({
 "tot_tke"
});

std::vector<std::string> profs_dycoms({
"00rtot", "rliq", "thl", "wvar", 
"w3rd", "prflux"
,"clfrac"
//, "N_c", 
,"cl_nc"
,"sat_RH"
,"rad_flx"
//, "nc_up" 
//,"sat_RH_up"
//, "act_conc_up" 
//, "nc_down" 
}); // rtot has to be first

std::vector<std::string> profs_sgs({
 "sgs_tke"
,"k_m"
,"sgs_tht_flux"
,"sgs_rv_flux"
//,"sgs_rc_flux"
,"sgs_u_flux"
});

std::vector<std::string> profs_moist_thermal({
}); // rtot has to be first


std::vector<std::string> fields_dycoms({
"rl", "nc",
 "rr", "nr",
"ef", "na", 
"th", "rv",     
"u", "w", 
"sd_conc",//, "r_dry", 
"RH", "supersat",
"lib_pres", "lib_temp"
});

std::vector<std::string> fields_moist_thermal({
//"mrk", "vel_div",
"vel_div",
"rl", "nc",
// "rr", "nr",
//"ef", "na", 
"th", "rv",     
"u", "w", 
"sd_conc",//, "r_dry", 
"RH", "supersat"
});

class Plots
{
  public:
    std::vector<std::string> series;
    std::vector<std::string> profs;
    std::vector<std::string> fields;

  Plots(const std::string &type, bool sgs):
    series(type == "dycoms" ? series_dycoms : series_moist_thermal),
    profs(type == "dycoms" ? profs_dycoms : profs_moist_thermal),
    fields(type == "dycoms" ? fields_dycoms : fields_moist_thermal)
  {
    if (sgs)
    {
      profs.insert(profs.end(), profs_sgs.begin(), profs_sgs.end());
      series.insert(series.end(), series_sgs.begin(), series_sgs.end());
    }
  }
};
