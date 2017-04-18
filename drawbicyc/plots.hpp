#pragma once

const std::set<std::string> series_dycoms({
"wvarmax", "clfrac", "lwp", "er",
 "surf_precip", "mass_dry", "acc_precip",
 "cl_nc",
 "sd_conc_avg", "sd_conc_std_dev",
 "tot_water"
});

const std::set<std::string> series_moist_thermal({
 "mass_dry",
"ract_com",
 "ract_avg", "ract_std_dev",
"com_mom0","com_mom1","com_mom2", // higher moments need lower ones enabled!!
 "sd_conc_avg", "sd_conc_std_dev",
"RH_max"
"th_com"
});

std::set<std::string> profs_dycoms({
"00rtot", "rliq", "thl", "wvar", 
"w3rd", "prflux", "act_conc", 
"clfrac", "N_c", 
"sat_RH"
}); // rtot has to be first

std::set<std::string> profs_moist_thermal({
}); // rtot has to be first


std::set<std::string> fields_dycoms({
"rl", "nc",
 "rr", "nr",
"ef", "na", 
"th", "rv",     
"u", "w", 
"sd_conc",//, "r_dry", 
"RH", "supersat"
});

std::set<std::string> fields_moist_thermal({
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
    const std::set<std::string> series;
    const std::set<std::string> profs;
    const std::set<std::string> fields;

  Plots(const std::string &type):
    series(type == "dycoms" ? series_dycoms : series_moist_thermal),
    profs(type == "dycoms" ? profs_dycoms : profs_moist_thermal),
    fields(type == "dycoms" ? fields_dycoms : fields_moist_thermal)
  {}
};
