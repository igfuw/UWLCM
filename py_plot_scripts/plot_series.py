import h5py
import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys

from plot_hlpr import read_data

###########
# plot time series of variables:
# LWP, cloud fraction, cloud thickness
###########

def main():
    path = sys.argv[1]
    test = sys.argv[2]
    folder = "../"+path+"/"+test+"/"
    
    calc_stats(folder)

    plot_lwp(folder)
    plot_cc(folder)
    plot_cthickness(folder)
    plot_precip(folder)

    # folders = ['../DYCOMS/lgr_r0.1_n100/', '../DYCOMS/lgr_r0.01_n100/', '../DYCOMS/lgr_r0.1_n1000/', '../DYCOMS/lgr_r0.01_n1000/',]
    # plot_mult_lwp(folders)
    # plot_mult_cc(folders)
    # plot_mult_cthickness(folders)

############
# functions
############

def calc_stats(folder):
    filelist = glob.glob(folder+"timestep*"+"*.h5")
    list.sort(filelist, key=lambda x: int(x.split("step")[-1].split(".h5")[0]))

    # load const
    with h5py.File(folder+"const.h5","r") as f:
        T = np.array(f["T"])
        Z = np.array(f["Y"])
        rhod = np.array(f["rhod"])
    
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

    if not os.path.exists(folder+'plots/'):
        os.makedirs(folder+'plots/')
    np.savetxt(folder+"plots/dat.csv", dat, delimiter="\t",fmt="%.2e",header=
        "Liquid water path\nTime, Mean LWP (g/m^2), Std LWP (g/m^2), Mean Cloud Cover, Std Cloud Cover, Mean Thickness (m), Mean Surface Precip (mm/d)")

def calc_lwp(file, rhod, z):
    cwc = read_data(file,'cwc')  # (g/kg)
    rwc = read_data(file,'rwc')  # (g/kg)
    rv = read_data(file,'rv')   # (g/kg)

    rho = rhod * (1 + rv/1e3 + cwc/1e3 + rwc/1e3)

    dz = np.diff(z)[0:-1,:]
    integrand = rho * cwc * dz
    lwp = np.sum(integrand, axis = 1)
    return lwp

def calc_ccover(file, rhod):
    nc = read_data(file,'nc')
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

    nc = read_data(file,'nc')
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

def plot_lwp(folder):
    path = folder+"plots/dat.csv"
    dat = np.genfromtxt(path,delimiter="\t")
    time = dat[1:,0]
    mean = dat[1:,1]
    std = dat[1:,2]

    plt.plot(time/3600., mean, "bo-")
    plt.fill_between(time/3600., mean-std, mean+std, facecolor="b", alpha=0.5)

    # set plot attributes
    fs = 12
    plt.xlabel("Time (hours)",fontsize=fs)
    plt.ylabel("LWP (g/m$^2$)",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title("Liquid water path", fontsize=fs)

    plt.tight_layout()
    plt.savefig(folder+"plots/lwp.png",dpi=300)
    plt.close()

def plot_cc(folder):
    path = folder+"plots/dat.csv"
    dat = np.genfromtxt(path,delimiter="\t")
    time = dat[1:,0]
    mean = dat[1:,3]
    std = dat[1:,4]

    plt.plot(time/3600., mean, "bo-")
    plt.fill_between(time/3600., mean-std, mean+std, facecolor="b", alpha=0.5)

    # set plot attributes
    fs = 12
    plt.xlabel("Time (hours)",fontsize=fs)
    plt.ylabel("Cloud cover",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title("Cloud cover", fontsize=fs)

    plt.tight_layout()
    plt.savefig(folder+"plots/ccover.png",dpi=300)
    plt.close()

def plot_cthickness(folder):
    path = folder+"plots/dat.csv"
    dat = np.genfromtxt(path,delimiter="\t")
    time = dat[1:,0]
    cthick = dat[1:,5]

    plt.plot(time/3600., cthick, "bo-")

    # set plot attributes
    fs = 12
    plt.xlabel("Time (hours)",fontsize=fs)
    plt.ylabel("Thickness (m)",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title("Cloud thickness", fontsize=fs)

    plt.tight_layout()
    plt.savefig(folder+"plots/cthick.png",dpi=300)
    plt.close()

def plot_precip(folder):
    path = folder+"plots/dat.csv"
    dat = np.genfromtxt(path,delimiter="\t")
    time = dat[1:,0]
    precip = dat[1:,6]

    plt.plot(time/3600., precip, "bo-")

    # set plot attributes
    fs = 12
    plt.xlabel("Time (hours)",fontsize=fs)
    plt.ylabel("Surface precip. (mm day$^{-1}$)",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title("Surface precipitation", fontsize=fs)

    plt.tight_layout()
    plt.savefig(folder+"plots/surf_precip.png",dpi=300)
    plt.close()

def plot_mult_lwp(folders):
    for folder in folders:
        path = folder+"plots/dat.csv"
        dat = np.genfromtxt(path,delimiter="\t")
        time = dat[1:,0]
        mean = dat[1:,1]
        #std = dat[1:,2]

        name = folder.split("/")[2]
        plt.plot(time/3600., mean, "o-",label=name)
        #plt.fill_between(time/3600., mean-std, mean+std, facecolor="b", alpha=0.5)

    # set plot attributes
    fs = 12
    plt.xlabel("Time (hours)",fontsize=fs)
    plt.ylabel("LWP (g/m$^2$)",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title("Liquid water path", fontsize=fs)
    plt.legend(loc=4)

    plt.tight_layout()
    plt.savefig("mult_lwp.png",dpi=300)
    plt.close()

def plot_mult_cc(folders):
    for folder in folders:
        path = folder+"plots/dat.csv"
        dat = np.genfromtxt(path,delimiter="\t")
        time = dat[1:,0]
        mean = dat[1:,3]
        std = dat[1:,4]

        name = folder.split("/")[2]
        plt.plot(time/3600., mean, "o-",label=name)
        #plt.fill_between(time/3600., mean-std, mean+std, facecolor="b", alpha=0.5)

    # set plot attributes
    fs = 12
    plt.xlabel("Time (hours)",fontsize=fs)
    plt.ylabel("Cloud cover (%)",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title("Cloud fraction", fontsize=fs)
    plt.legend(loc=3)

    plt.tight_layout()
    plt.savefig("mult_cf.png",dpi=300)
    plt.close()

def plot_mult_cthickness(folders):
    for folder in folders:
        path = folder+"plots/dat.csv"
        dat = np.genfromtxt(path,delimiter="\t")
        time = dat[1:,0]
        cthick = dat[1:,5]

        name = folder.split("/")[2]
        plt.plot(time/3600., cthick, "o-", label=name)

    # set plot attributes
    fs = 12
    plt.xlabel("Time (hours)",fontsize=fs)
    plt.ylabel("Thickness (m)",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title("Cloud thickness", fontsize=fs)
    plt.legend(loc=4)

    plt.tight_layout()
    plt.savefig("mult_cthick.png",dpi=300)
    plt.close()

if __name__ == "__main__":
    main()