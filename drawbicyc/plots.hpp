#pragma once

const std::vector<std::string> series_dycoms({
"wvarmax", "clfrac", "lwp", "er",
 "surf_precip", 
//"mass_dry", 
 "acc_precip",
 "cl_nc",
 "cloud_base",
 "cl_gccn_conc", "gccn_conc", "cl_gccn_meanr"
,"cl_avg_cloud_rad"
// "sd_conc_avg", "sd_conc_std_dev",
// "tot_water"
});

const std::vector<std::string> series_rico({
 "clfrac", "lwp",
 "surf_precip", 
 "acc_precip",
 "cl_nc",
 "cloud_base",
 "surf_flux_latent",
 "surf_flux_sensible"
//"mass_dry", 
// "cl_gccn_conc", "gccn_conc", "cl_gccn_meanr"
//,"cl_avg_cloud_rad"
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

std::vector<std::string> profs_rico({
"00rtot", "rliq", "thl", "wvar", 
 "prflux"
,"clfrac"
//, "N_c", 
,"cl_nc"
,"u", "v",
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
//"rl", "nc",
// "rr", "nr",
//"ef", "na", 
//"th", "rv",     
//"u", "w", 
//"sd_conc",//, "r_dry", 
//"RH", "supersat",
//"lib_pres", "lib_temp"
"gccn_conc",
"gccn_mean_rw"
});

std::vector<std::string> fields_rico({
"rl", "nc",
 "rr", "nr",
//"ef", "na", 
"th", "rv",     
"u", "w", 
//"sd_conc",//, "r_dry", 
//"RH", "supersat",
//"lib_pres", "lib_temp"
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

  Plots(const std::string &type, bool sgs)
  {
    if(type == "dycoms") { 
      profs.insert(profs.end(), profs_dycoms.begin(), profs_dycoms.end());
      series.insert(series.end(), series_dycoms.begin(), series_dycoms.end());
      fields.insert(fields.end(), fields_dycoms.begin(), fields_dycoms.end());
    }
    else if(type == "rico") { 
      profs.insert(profs.end(), profs_rico.begin(), profs_rico.end());
      series.insert(series.end(), series_rico.begin(), series_rico.end());
      fields.insert(fields.end(), fields_rico.begin(), fields_rico.end());
    }
    else if(type == "moist_thermal") { 
      profs.insert(profs.end(), profs_moist_thermal.begin(), profs_moist_thermal.end());
      series.insert(series.end(), series_moist_thermal.begin(), series_moist_thermal.end());
      fields.insert(fields.end(), fields_moist_thermal.begin(), fields_moist_thermal.end());
    }
    else
      throw std::runtime_error("drawbicyc Plots.hpp: unknown 'type'.");
    
    if (sgs)
    {
      profs.insert(profs.end(), profs_sgs.begin(), profs_sgs.end());
      series.insert(series.end(), series_sgs.begin(), series_sgs.end());
    }
  }
};
