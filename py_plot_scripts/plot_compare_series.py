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

    plot_all_compare(path, exps)

def plot_all_compare(path, exps):
    folders = [path + x for x in exps]
    files = [x + "/ensemble_mean_dat.csv" for x in folders]

    fig = plt.figure(figsize=(18,12))
    for file in files:
        dat = np.genfromtxt(file,delimiter="\t")
        time = dat[1:,0]

        exp_name = file.split("/")[3]
        kappa = exp_name.split("_")[-1]
        name = "$\kappa = $" + kappa

        plt.subplot(221)
        mean = dat[1:,1]
        std = dat[1:,2]
        plt.plot(time/3600., mean, "o-",label=name)
        plt.fill_between(time/3600., mean-std, mean+std, alpha=0.5)

        plt.subplot(222)
        mean = dat[1:,3]
        std = dat[1:,4]
        plt.plot(time/3600., mean, "o-",label=name)
        plt.fill_between(time/3600., mean-std, mean+std, alpha=0.5)

        plt.subplot(223)
        mean = dat[1:,5]
        std = dat[1:,6]
        plt.plot(time/3600., mean, "o-",label=name)
        plt.fill_between(time/3600., mean-std, mean+std, alpha=0.5)
        
        plt.subplot(224)
        mean = dat[1:,7]
        std = dat[1:,8]
        plt.plot(time/3600., mean, "o-",label=name)
        plt.fill_between(time/3600., mean-std, mean+std, alpha=0.5)

    # set subplot attributes
    plt.subplot(221)
    plt.xlabel("Time (hours)")
    plt.ylabel("Liquid water path (g/m$^2$)")
    plt.legend(loc=3)

    plt.subplot(222)
    plt.xlabel("Time (hours)")
    plt.ylabel("Cloud cover")

    plt.subplot(223)
    plt.xlabel("Time (hours)")
    plt.ylabel("Cloud thickness (m)")

    plt.subplot(224)
    plt.xlabel("Time (hours)")
    plt.ylabel("Surface precipitation (mm/day)")

    # set supplot title
    plt.suptitle("mean_rd = 0.01e-6, n_stp = 100e6")

    fs = 20
    plt.rcParams.update({'font.size': fs})
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)

    if not os.path.exists(path+"comparison_plots/"):
        os.makedirs(path+"comparison_plots/")
    plt.savefig(path+"comparison_plots/var_kappa.png",dpi=300)
    plt.close()

if __name__ == "__main__":
    main()
