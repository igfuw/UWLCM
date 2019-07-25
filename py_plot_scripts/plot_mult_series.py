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
    path = '../../output/mu_0.01e-6_kappa_1.3/'
    folder_list = ['exp1/', 'exp2/', 'exp3/', 'exp4/', 'exp5/']

    calc_avgs(path, folder_list)

    plot_mult_lwp(path, folder_list)
    plot_mult_cc(path, folder_list)
    plot_mult_cthickness(path, folder_list)
    plot_mult_precip(path, folder_list)

def calc_avgs(path, folder_list):
    folders = [path + x for x in folder_list]
    file = folders[0] + "plots/dat.csv"
    dat = np.genfromtxt(file,delimiter="\t")
    time = dat[1:,0]
    
    lwp_arr = np.zeros((len(folder_list), len(time)))
    cc_arr = np.zeros((len(folder_list), len(time)))
    cthick_arr = np.zeros((len(folder_list), len(time)))
    precip_arr = np.zeros((len(folder_list), len(time)))

    for i,folder in enumerate(folders):
        file = folder+"plots/dat.csv"
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
    np.savetxt(path+"ensemble_mean_dat.csv", dat, delimiter="\t",fmt="%.2e",header="Time, LWP (g/m^2), std, Cloud Cover, std, Cloud Thickness (m), std, Surf Precip (mm/d), std")
    

def plot_mult_lwp(path, folder_list):
    folders = [path + x for x in folder_list]
    for folder in folders:
        file = folder+"plots/dat.csv"
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
        file = folder+"plots/dat.csv"
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
        file = folder+"plots/dat.csv"
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
        file = folder+"plots/dat.csv"
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
