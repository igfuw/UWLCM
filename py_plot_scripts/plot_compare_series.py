import h5py
import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys

def main():
    path="../../output/"
    exps = ["mu_0.01e-6_kappa_0.1","mu_0.01e-6_kappa_1.3"]

    plot_lwp_compare(path, exps)

def plot_lwp_compare(path, exps):
    folders = [path + x for x in exps]
    files = [x + "/ensemble_mean_dat.csv" for x in folders]
    for file in files:
        dat = np.genfromtxt(file,delimiter="\t")
        time = dat[1:,0]
        lwp = dat[1:,1]
        std = dat[1:,2]

        name = file.split("/")[3]
        plt.plot(time/3600., lwp, "o-",label=name)
        plt.fill_between(time/3600., lwp-std, lwp+std, alpha=0.5)

    # set plot attributes
    fs = 12
    plt.xlabel("Time (hours)",fontsize=fs)
    plt.ylabel("LWP (g/m$^2$)",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title("Liquid water path", fontsize=fs)
    plt.legend(loc=3)

    plt.tight_layout()

    if not os.path.exists(path+"comparison_plots/"):
        os.makedirs(path+"comparison_plots/")
    plt.savefig(path+"comparison_plots/lwp.png",dpi=300)
    plt.close()

if __name__ == "__main__":
    main()
