import h5py
import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import sys

from read_data import read_hdf

def main():
    path = sys.argv[1]
    test = sys.argv[2]
    folder = "../"+path+"/"+test+"/"

    calc_profs(folder)
    plot_cf_z(folder)

def calc_profs(folder):
    filelist = glob.glob(folder+"timestep*"+"*.h5")
    list.sort(filelist, key=lambda x: int(x.split("step")[-1].split(".h5")[0]))

    # load const
    with h5py.File(folder+"const.h5","r") as f:
        T = np.array(f["T"])
        Z = np.array(f["Y"])
        rhod = np.array(f["rhod"])
    
    Z = Z[0,:-1]
    cf_z = np.zeros((len(Z)+1,len(T)+1))
    cf_z[0,0] = np.nan
    cf_z[0,1:] = T
    cf_z[1:,0] = Z
    for t,file in enumerate(filelist):
        cf_z[1:,t+1] = calc_cf_z(file, rhod)

    if not os.path.exists(folder+'plots/'):
        os.makedirs(folder+'plots/')
    np.savetxt(folder+"plots/cf_z.csv", cf_z, delimiter="\t",fmt="%.2f",header="Z, CF @ Time T")

def calc_cf_z(file, rhod):
    nc = read_hdf(file,'nc')
    rhod = np.tile(rhod, np.shape(nc)[0]).reshape(np.shape(nc))
    nc *= rhod # (1/cm^3)
    cloudy_cutoff = 50 # cloudy if nc > X mg^-1
    cloudy = (nc > cloudy_cutoff).astype(int)
    cf_z = np.mean(cloudy, axis=0)
    return cf_z

def plot_cf_z(folder):
    path = folder+"plots/cf_z.csv"
    dat = np.genfromtxt(path,delimiter="\t")
    
    time = dat[0,1:]
    z = dat[1:,0]

    for i,t in enumerate(time):
        cf = dat[1:,i+1]
        plt.plot(cf,z,"o-",c=cm.plasma(t/np.max(time)),label="{:.2f}".format(t/3600.))

    # set plot attributes
    fs = 12
    plt.xlabel("Cloud fraction (%)",fontsize=fs)
    plt.ylabel("Altitude (m)",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title("Cloud fraction", fontsize=fs)
    plt.legend(loc="upper right",title="Time (hours)",ncol=3)

    plt.tight_layout()
    plt.savefig(folder+"plots/cf_z.png",dpi=300)
    plt.close()

if __name__ == "__main__":
    main()