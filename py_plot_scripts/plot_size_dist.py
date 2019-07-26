import h5py
import numpy as np
import glob
from matplotlib import cm
import matplotlib.pyplot as plt
import os
import sys

from read_data import read_bins

def main():
    path = sys.argv[1]
    test = sys.argv[2]
    folder = "../"+path+"/"+test+"/"

    calc_sizes(folder)
    plot_wet_dist(folder)
    plot_wet_dry_dist(folder)

def calc_sizes(folder):
    filelist = glob.glob(folder+"timestep*"+"*.h5")
    list.sort(filelist, key=lambda x: int(x.split("step")[-1].split(".h5")[0]))

    # load const
    with h5py.File(folder+"const.h5","r") as f:
        T = np.array(f["T"])
    
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
    
    if not os.path.exists(folder+'plots/'):
        os.makedirs(folder+'plots/')
    np.savetxt(folder+"plots/wet_dist.csv", dat, delimiter="\t",fmt="%.3e",
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
    
    np.savetxt(folder+"plots/dry_dist.csv", dat, delimiter="\t",fmt="%.3e",
        header="Radius (um), Bin width (um), Num. conc. (mg-1) @ Time T (hours)")

def calc_dist(file, r, bin_type):
    n = np.zeros(len(r))
    with h5py.File(file,"r") as f:
        for bin_num in np.arange(len(n)):
            var = "r{}_rng{:03d}_mom0".format(bin_type, bin_num) # conc = number per mg of air
            dat = 1e-6 * np.array(f[var])
            avg = np.mean(dat)
            n[bin_num] = avg
    return n

def plot_wet_dist(folder):
    path = folder+"plots/wet_dist.csv"
    dat = np.genfromtxt(path,delimiter="\t")
    
    time = dat[0,2:]
    r = dat[1:,0]
    dr = dat[1:,1]
    
    for i,t in enumerate(time):
        n = dat[1:,i+2]
        plt.loglog(r,n/dr,"o-",c=cm.plasma(t/np.max(time)),label="{:.2f}".format(t))

        # set plot attributes
    fs = 12
    plt.xlabel(r"Wet radius ($\mu$m)",fontsize=fs)
    plt.ylabel(r"dn/dr (mg$^{-1}$ $\mu$m$^{-1}$)",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.xlim([1e-3,1e2])
    plt.ylim([1e-2,1e5])
    plt.title("Droplet size distribution", fontsize=fs)
    plt.legend(loc="upper right",title="Time (hours)",ncol=3)

    plt.tight_layout()
    plt.savefig(folder+"plots/wet_dist.png",dpi=300)
    plt.close()

def plot_wet_dry_dist(folder):
    # plot dry
    path = folder+"plots/dry_dist.csv"
    dat = np.genfromtxt(path,delimiter="\t")
    r = dat[1:,0]
    dr = dat[1:,1]
    n = dat[1:,-1]
    plt.loglog(r,n/dr,"ks--",label="dry")

    # plot wet
    path = folder+"plots/wet_dist.csv"
    dat = np.genfromtxt(path,delimiter="\t")
    r = dat[1:,0]
    dr = dat[1:,1]
    n = dat[1:,-1]
    plt.loglog(r,n/dr,"bo-",label="wet")

    # set plot attributes
    fs = 12
    plt.xlabel(r"Particle radius ($\mu$m)",fontsize=fs)
    plt.ylabel(r"dn/dr (mg$^{-1}$ $\mu$m$^{-1}$)",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.xlim([1e-3,1e2])
    plt.ylim([1e-2,1e5])
    plt.title("Droplet size distribution", fontsize=fs)
    plt.legend(loc="upper right")

    # save plot
    plt.tight_layout()
    plt.savefig(folder+"plots/wet_dry_dist.png",dpi=300)
    plt.close()

if __name__ == "__main__":
    main()