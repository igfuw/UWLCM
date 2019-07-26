import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

def plot_mult_series_all(path, folder_list):
    plt.figure(figsize=(18,12))

    folders = [path + x for x in folder_list]

    for i,folder in enumerate(folders):
        name = folder_list[i]
        file = folder+"plots/dat_series.csv"
        dat = np.genfromtxt(file,delimiter="\t")
        time = dat[1:,0]
        lwp = dat[1:,1]
        lwp_std = dat[1:,2]
        cc = dat[1:,3]
        cc_std = dat[1:,4]
        cthick = dat[1:,5]
        precip = dat[1:,6]

        plt.subplot(221)
        plt.plot(time/3600., lwp, "o-", label=name)
        plt.fill_between(time/3600., lwp-lwp_std, lwp+lwp_std, alpha=0.5)

        plt.subplot(222)
        plt.plot(time/3600., cc, "o-", label=name)
        plt.fill_between(time/3600., cc-cc_std, cc+cc_std, alpha=0.5)

        plt.subplot(223)
        plt.plot(time/3600., cthick, "o-", label=name)
        
        plt.subplot(224)
        plt.plot(time/3600., precip, "o-", label=name)

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
    plt.ylabel("Accumulated surface precipitation (mm/day)")

    # set supplot title
    plt.suptitle(path)
    fs = 20
    plt.rcParams.update({'font.size': fs})
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.show()

    plt.savefig(path+"plots/mult_series_all.png",dpi=300)
    plt.close()

def plot_series_all(folder):
    plt.figure(figsize=(18,12))

    file = folder+"plots/dat_series.csv"
    dat = np.genfromtxt(file,delimiter="\t")
    time = dat[1:,0]
    lwp = dat[1:,1]
    lwp_std = dat[1:,2]
    cc = dat[1:,3]
    cc_std = dat[1:,4]
    cthick = dat[1:,5]
    precip = dat[1:,6]

    plt.subplot(221)
    plt.plot(time/3600., lwp, "o-")
    plt.fill_between(time/3600., lwp-lwp_std, lwp+lwp_std, alpha=0.5)

    plt.subplot(222)
    plt.plot(time/3600., cc, "o-")
    plt.fill_between(time/3600., cc-cc_std, cc+cc_std, alpha=0.5)

    plt.subplot(223)
    plt.plot(time/3600., cthick, "o-")
    
    plt.subplot(224)
    plt.plot(time/3600., precip, "o-")

    # set subplot attributes
    plt.subplot(221)
    plt.xlabel("Time (hours)")
    plt.ylabel("Liquid water path (g/m$^2$)")

    plt.subplot(222)
    plt.xlabel("Time (hours)")
    plt.ylabel("Cloud cover")

    plt.subplot(223)
    plt.xlabel("Time (hours)")
    plt.ylabel("Cloud thickness (m)")

    plt.subplot(224)
    plt.xlabel("Time (hours)")
    plt.ylabel("Accumulated surface precipitation (mm/day)")

    # set supplot title
    plt.suptitle(folder)
    fs = 20
    plt.rcParams.update({'font.size': fs})
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.show()

    plt.savefig(folder+"plots/series_all.png",dpi=300)
    plt.close()

def plot_cf_z(folder):
    path = folder+"plots/dat_cf_prof.csv"
    dat = np.genfromtxt(path,delimiter="\t")
    
    time = dat[0,1:]
    z = dat[1:,0]

    for i,t in enumerate(time):
        if np.mod(t,3600) == 0:
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
    plt.show()

    # save fig
    plt.savefig(folder+"plots/cf_z.png",dpi=300)
    plt.close()

def plot_wet_dist_evo(folder):
    path = folder+"plots/wet_dist.csv"
    dat = np.genfromtxt(path,delimiter="\t")
    
    time = dat[0,2:]
    r = dat[1:,0]
    dr = dat[1:,1]
    
    for i,t in enumerate(time):
        if np.mod(t,3600) == 0:
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
    plt.show()

    # save fig
    plt.savefig(folder+"plots/wet_dist_evo.png",dpi=300)
    plt.close()

def plot_wet_dry_dist(folder):
    # plot init dry
    path = folder+"plots/dat_dry_dist.csv"
    dat = np.genfromtxt(path,delimiter="\t")
    r = dat[1:,0]
    dr = dat[1:,1]
    n = dat[1:,2]
    plt.loglog(r,n/dr,"ks--",label="init dry")

    # plot final dry
    n = dat[1:,-1]
    plt.loglog(r,n/dr,"ks-",label="final dry")

    # plot init wet
    path = folder+"plots/dat_wet_dist.csv"
    dat = np.genfromtxt(path,delimiter="\t")
    r = dat[1:,0]
    dr = dat[1:,1]
    n = dat[1:,2]
    plt.loglog(r,n/dr,"bo--",label="init wet")

    # plot final wet
    n = dat[1:,-1]
    plt.loglog(r,n/dr,"bo-",label="final wet")

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
    plt.tight_layout()
    plt.show()

    # save plot
    plt.savefig(folder+"plots/wet_dry_dist.png",dpi=300)
    plt.close()