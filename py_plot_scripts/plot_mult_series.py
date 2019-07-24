import h5py
import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys

###########
# plot time series of variables:
# LWP, cloud fraction, cloud thickness, surface precip
###########

def main():
    path = '../output_lgr/dycoms/6hrs_highres/'
    folder_list = ['exp1/', 'exp2/', 'exp3/', 'exp4/', 'exp5/']
    plot_mult_lwp(path, folder_list)
    plot_mult_cc(path, folder_list)
    plot_mult_cthickness(path, folder_list)
    plot_mult_precip(path, folder_list)

def plot_mult_lwp(path, folder_list):
    folders = [path + x for x in folder_list]
    for folder in folders:
        file = folder+"dat.csv"
        dat = np.genfromtxt(file,delimiter="\t")
        time = dat[1:,0]
        mean = dat[1:,1]

        name = folder.split("/")[4]
        plt.plot(time/3600., mean, "o-",label=name)

    # set plot attributes
    fs = 12
    plt.xlabel("Time (hours)",fontsize=fs)
    plt.ylabel("LWP (g/m$^2$)",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title("Liquid water path", fontsize=fs)
    plt.legend(loc=4)

    plt.tight_layout()
    plt.savefig(path+"mult_lwp.png",dpi=300)
    plt.close()
    
def plot_mult_cc(path, folder_list):
    folders = [path + x for x in folder_list]
    for folder in folders:
        file = folder+"dat.csv"
        dat = np.genfromtxt(file,delimiter="\t")
        time = dat[1:,0]
        mean = dat[1:,3]
        #std = dat[1:,4]

        name = folder.split("/")[4]
        plt.plot(time/3600., mean, "o-",label=name)
        #plt.fill_between(time/3600., mean-std, mean+std, facecolor="b", alpha=0.5)

    # set plot attributes
    fs = 12
    plt.xlabel("Time (hours)",fontsize=fs)
    plt.ylabel("Cloud cover (%)",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title("Cloud cover", fontsize=fs)
    plt.legend(loc=3)

    plt.tight_layout()
    plt.savefig(path+"mult_cc.png",dpi=300)
    plt.close()

def plot_mult_cthickness(path, folder_list):
    folders = [path + x for x in folder_list]
    for folder in folders:
        file = folder+"dat.csv"
        dat = np.genfromtxt(file,delimiter="\t")
        time = dat[1:,0]
        cthick = dat[1:,5]

        name = folder.split("/")[4]
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
    plt.savefig(path+"mult_cthick.png",dpi=300)
    plt.close()

def plot_mult_precip(path, folder_list):
    folders = [path + x for x in folder_list]
    for folder in folders:
        file = folder+"dat.csv"
        dat = np.genfromtxt(file,delimiter="\t")
        time = dat[1:,0]
        precip = dat[1:,6]

        name = folder.split("/")[4]
        plt.plot(time/3600., precip, "o-", label=name)

    # set plot attributes
    fs = 12
    plt.xlabel("Time (hours)",fontsize=fs)
    plt.ylabel("Surface precip. (mm/day)",fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title("Surface precipitation", fontsize=fs)
    plt.legend(loc=4)

    plt.tight_layout()
    plt.savefig(path+"mult_precip.png",dpi=300)
    plt.close()

if __name__ == "__main__":
    main()