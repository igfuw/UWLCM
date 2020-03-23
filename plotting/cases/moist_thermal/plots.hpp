#pragma once

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


std::vector<std::string> profs_moist_thermal({
}); // rtot has to be first


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
