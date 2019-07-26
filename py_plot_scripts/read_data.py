import h5py
import numpy as np

def read_hdf(file,var):
    water_density = 1e3 # (kg/m^3)
    # load data
    with h5py.File(file,'r') as f:
        if var == 'cwc':
            dat = np.array(f['cloud_rw_mom3']) * 4/3 * np.pi * water_density * 1e3 # cloud water content (g/kg)
        elif var == 'rwc':
            dat = np.array(f['rain_rw_mom3']) * 4/3 * np.pi * water_density * 1e3 # rain water content (g/kg)
        elif var == 'nc':
            dat = np.array(f['cloud_rw_mom0']) * 1e-6 # cloud particle concentration (1/mg)
        elif var == 'nr':
            dat = np.array(f['rain_rw_mom0']) * 1e-6 # rain particle concentration (1/mg)
        elif var == 'rv':
            dat = np.array(f['rv']) * 1000. # water vapor content (g/kg)
        elif var == 'th':
            dat = np.array(f['th'])     # potential temperature (K)
        elif var == 'u':
            dat = np.array(f['u'])      # horizontal velocity (m/s)
        elif var == 'w':
            dat = np.array(f['w'])      # vertical velocity (m/s)
        elif var == 'precip_rate':
            dat = np.array(f['precip_rate']) * 4/3 * np.pi * water_density   # precip rate (m/s)
        # elif var == 'ef':
        #     dat = np.array(f['rw_rng000_mom3']) / np.array(f['rw_rng000_mom2']) * 1e6 # effective radius (um)
        else:
            dat = np.array(f[var])
    return dat

def read_bins(folder, bin_type):
    with open(folder+"bin_strs.txt","r") as f:
        lines = f.read().splitlines()
    lines = np.array(lines)

    line_num = np.where(lines == bin_type+"_bins_str:")[0][0] + 1
    bins_str = lines[line_num] 
    
    x = bins_str.replace(';',':')
    x = x.split(":")
    x = x[::2]
    left_edges = np.array([float(i) for i in x])
    centers = np.sqrt(left_edges[0:-1]*left_edges[1:])
    widths = left_edges[1:] - left_edges[0:-1]
    return (centers, widths)