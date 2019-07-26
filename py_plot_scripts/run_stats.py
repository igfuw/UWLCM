import numpy as np
import glob
import os
import sys

from read_data import read_hdf, read_bins
from stats_hlpr import calc_dist, calc_lwp, calc_ccover, calc_cthickness, calc_surf_precip, calc_cf_z

def main():
    path = sys.argv[1]
    test = sys.argv[2]
    folder = path+"/"+test+"/"
    print(folder)

    if not os.path.exists(folder+'plots/'):
        os.makedirs(folder+'plots/')

    calc_sizes(folder)
    calc_series(folder)
    calc_profs(folder)

    folder_list = ['exp1/', 'exp2/', 'exp3/', 'exp4/', 'exp5/']
    calc_series_avgs(path, folder_list)

def calc_sizes(folder):
    filelist = glob.glob(folder+"timestep*"+"*.h5")
    list.sort(filelist, key=lambda x: int(x.split("step")[-1].split(".h5")[0]))

    # load const
    T = read_hdf(folder+"const.h5","T")
    
    # wet
    (r, dr) = read_bins(folder, "wet")
    dat = np.zeros((len(r)+1,len(T)+2))
    dat[0,0] = np.nan
    dat[0,1] = np.nan
    dat[1:,0] = r / 1e-6
    dat[1:,1] = dr / 1e-6
    dat[0,2:] = T / 3600.
    for t,file in enumerate(filelist):
        n = calc_dist(file, r, "w")
        dat[1:,t+2] = n
    
    np.savetxt(folder+"plots/dat_wet_dist.csv", dat, delimiter="\t",fmt="%.3e",
        header="Radius (um), Bin width (um), Num. conc. (mg-1) @ Time T (hours)")

    # dry
    (r, dr) = read_bins(folder, "dry")
    dat = np.zeros((len(r)+1,len(T)+2))
    dat[0,0] = np.nan
    dat[0,1] = np.nan
    dat[1:,0] = r / 1e-6
    dat[1:,1] = dr / 1e-6
    dat[0,2:] = T / 3600.
    for t,file in enumerate(filelist):
        n = calc_dist(file, r, "d")
        dat[1:,t+2] = n
    
    np.savetxt(folder+"plots/dat_dry_dist.csv", dat, delimiter="\t",fmt="%.3e",
        header="Radius (um), Bin width (um), Num. conc. (mg-1) @ Time T (hours)")

def calc_series(folder):
    filelist = glob.glob(folder+"timestep*"+"*.h5")
    list.sort(filelist, key=lambda x: int(x.split("step")[-1].split(".h5")[0]))

    # load const
    T = read_hdf(folder+"const.h5","T")
    Z = read_hdf(folder+"const.h5","Y")
    rhod = read_hdf(folder+"const.h5","rhod")

    # with h5py.File(folder+"const.h5","r") as f:
    #     T = np.array(f["T"])
    #     Z = np.array(f["Y"])
    #     rhod = np.array(f["rhod"])
    
    m_lwp = np.zeros(len(T))
    s_lwp = np.zeros(len(T))
    m_cc = np.zeros(len(T))
    s_cc = np.zeros(len(T))
    cthick = np.zeros(len(T))
    for t,file in enumerate(filelist):
        lwp = calc_lwp(file, rhod, Z)
        m_lwp[t] = np.mean(lwp)
        s_lwp[t] = np.std(lwp)
        cc = calc_ccover(file, rhod)
        m_cc[t] = np.mean(cc)
        s_cc[t] = np.std(cc)
        cthick[t] = calc_cthickness(file, rhod, Z)
    
    surf_precip = calc_surf_precip(folder,T)
    
    dat = np.array([T,m_lwp,s_lwp,m_cc,s_cc,cthick,surf_precip]).T

    np.savetxt(folder+"plots/dat_series.csv", dat, delimiter="\t",fmt="%.2e",header=
        "Liquid water path\nTime, Mean LWP (g/m^2), std, Mean Cloud Cover, std, Mean Thickness (m), Mean Surface Precip (mm/d)")

def calc_profs(folder):
    filelist = glob.glob(folder+"timestep*"+"*.h5")
    list.sort(filelist, key=lambda x: int(x.split("step")[-1].split(".h5")[0]))

    # load const
    T = read_hdf(folder+"const.h5","T")
    Z = read_hdf(folder+"const.h5","Y")
    rhod = read_hdf(folder+"const.h5","rhod")
    
    # with h5py.File(folder+"const.h5","r") as f:
    #     T = np.array(f["T"])
    #     Z = np.array(f["Y"])
    #     rhod = np.array(f["rhod"])
    
    Z = Z[0,:-1]
    dat = np.zeros((len(Z)+1,len(T)+1))
    dat[0,0] = np.nan
    dat[0,1:] = T
    dat[1:,0] = Z
    for t,file in enumerate(filelist):
        dat[1:,t+1] = calc_cf_z(file, rhod)

    np.savetxt(folder+"plots/dat_cf_prof.csv", dat, delimiter="\t",fmt="%.2f",header="Z, CF @ Time T")

def calc_series_avgs(path, folder_list):
    folders = [path + x for x in folder_list]
    file = folders[0] + "plots/dat_series.csv"
    dat = np.genfromtxt(file,delimiter="\t")
    time = dat[1:,0]
    
    lwp_arr = np.zeros((len(folder_list), len(time)))
    cc_arr = np.zeros((len(folder_list), len(time)))
    cthick_arr = np.zeros((len(folder_list), len(time)))
    precip_arr = np.zeros((len(folder_list), len(time)))

    for i,folder in enumerate(folders):
        file = folder+"plots/dat_series.csv"
        dat = np.genfromtxt(file,delimiter="\t")
        lwp_arr[i,:] = dat[1:,1]
        cc_arr[i,:] = dat[1:,3]
        cthick_arr[i,:] = dat[1:,5]
        precip_arr[i,:] = dat[1:,6]

    lwp = np.mean(lwp_arr,axis=0)
    lwp_std = np.std(lwp_arr,axis=0)
    cc = np.mean(cc_arr,axis=0)
    cc_std = np.std(cc_arr,axis=0)
    cthick = np.mean(cthick_arr,axis=0)
    cthick_std = np.std(cthick_arr,axis=0)
    precip = np.mean(precip_arr,axis=0)
    precip_std = np.std(precip_arr,axis=0)

    dat = np.array([time, lwp, lwp_std, cc, cc_std, cthick, cthick_std, precip, precip_std]).T
    np.savetxt(path+"ensemble_mean_dat.csv", dat, delimiter="\t",fmt="%.2e",header=
        "Time, LWP (g/m^2), std, Cloud Cover, std, Cloud Thickness (m), std, Surf Precip (mm/d), std")

if __name__ == "__main__":
    main()