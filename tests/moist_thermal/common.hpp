// Moist thermal comparison in setup from Clark Grabowski 1999 JAS
// only condensation and evaporation
// tested both for blk_1m and lgrngn microphysics
// second one takes long...

#include <unordered_map>
#include <string>

using std::unordered_map;
using std::string;

const unordered_map<string, string> outdir = { {"blk_1m", "tmp_out_blk_1m"}, {"lgrngn", "tmp_out_lgrngn"}};

const string opts_common = 
  "--outfreq=60 --nt=600 --spinup=0 --dt=1 --nx=181 --nz=121 --case=moist_thermal";
const unordered_map<string, string> opts_micro({
  {"blk_1m", "--micro=blk_1m --outdir="+outdir.at("blk_1m")+" --cond=true --cevp=true --revp=false --conv=false --accr=false --sedi=false"},
  {"lgrngn", "--micro=lgrngn --outdir="+outdir.at("lgrngn")+" --cond=true --adve=true --sedi=false --coal=false --backend=OpenMP --sd_conc=16 --rng_seed=44"}
});
