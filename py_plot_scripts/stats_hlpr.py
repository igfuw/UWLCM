import numpy as np
from read_data import read_hdf

def calc_dist(file, r, bin_type):
    n = np.zeros(len(r))
    for bin_num in np.arange(len(n)):
        varname = "r{}_rng{:03d}_mom0".format(bin_type, bin_num) # conc = number per mg of air
        dat = read_hdf(file,varname) * 1e-6
        avg = np.mean(dat)
        n[bin_num] = avg

    # with h5py.File(file,"r") as f:
    #     for bin_num in np.arange(len(n)):
    #         var = "r{}_rng{:03d}_mom0".format(bin_type, bin_num) # conc = number per mg of air
    #         dat = 1e-6 * np.array(f[var])
    #         avg = np.mean(dat)
    #         n[bin_num] = avg
    return n

def calc_lwp(file, rhod, z):
    cwc = read_hdf(file,'cwc')  # (g/kg)
    rwc = read_hdf(file,'rwc')  # (g/kg)
    rv = read_hdf(file,'rv')   # (g/kg)

    rho = rhod * (1 + rv/1e3 + cwc/1e3 + rwc/1e3)

    dz = np.diff(z)[0:-1,:]
    integrand = rho * cwc * dz
    lwp = np.sum(integrand, axis = 1)
    return lwp

def calc_ccover(file, rhod):
    nc = read_hdf(file,'nc')
    rhod = np.tile(rhod, np.shape(nc)[0]).reshape(np.shape(nc))
    nc *= rhod # (1/cm^3)
    cloudy_cutoff = 50 # cloudy if nc > X mg^-1
    cloudy = (nc > cloudy_cutoff).astype(int)
    col_cloudy = np.sum(cloudy, axis=1)
    col_cloudy = (col_cloudy > 0).astype(int)
    #cloud_cover = np.mean(col_cloudy)
    return col_cloudy

def calc_cthickness(file, rhod, Z):
    z = Z[0,:]

    nc = read_hdf(file,'nc')
    rhod = np.tile(rhod, np.shape(nc)[0]).reshape(np.shape(nc))
    nc *= rhod # (1/cm^3)
    cloudy_cutoff = 50 # cloudy if nc > cloudy_cutoff mg^-1
    cloudy = (nc > cloudy_cutoff).astype(int)
    lev_cf = np.mean(cloudy, axis=0)
    max_cf = np.max(lev_cf)
    if max_cf > 0:
        cf_cutoff = 0.2*max_cf # cloud if cf > cf_cutoff
        cloud_i = np.where(lev_cf > cf_cutoff)
        base_i = np.min(cloud_i)
        top_i = np.max(cloud_i)
        cthick = z[top_i] - z[base_i]
    else:
        cthick = 0
    return cthick

def calc_surf_precip(folder,T):
    with open(folder+"prec_vol.dat","r") as f:
        lines = f.read().splitlines()
    lines = np.array(lines)
    eighth = lines[8::11]

    avg_surf_precip = np.zeros(len(T))
    for t,line in enumerate(eighth):
        ind,precip = line.split(' ')
        avg_surf_precip[t] = precip

    return avg_surf_precip

def calc_cf_z(file, rhod):
    nc = read_hdf(file,'nc')
    rhod = np.tile(rhod, np.shape(nc)[0]).reshape(np.shape(nc))
    nc *= rhod # (1/cm^3)
    cloudy_cutoff = 50 # cloudy if nc > X mg^-1
    cloudy = (nc > cloudy_cutoff).astype(int)
    cf_z = np.mean(cloudy, axis=0)
    return cf_z