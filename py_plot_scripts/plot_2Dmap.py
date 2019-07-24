import h5py
import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import moviepy.editor as mpy
import sys

from plot_hlpr import read_data

###########
# plot 2d heat maps of variables:
# cloud/rain/vapour water content, cloud/rain particle concentration, effective radius
###########

def main():
    path = sys.argv[1]
    test = sys.argv[2]
    folder = "../"+path+'/'+test+'/'
    filelist = glob.glob(folder+'timestep*'+'*.h5')
    list.sort(filelist, key=lambda x: int(x.split('step')[-1].split('.h5')[0]))

    # load const
    with h5py.File(folder+'const.h5','r') as f:
        T = np.array(f['T'])
        X = np.array(f['X'])
        Z = np.array(f['Y'])

    vars = ['cwc','rwc','nc','nr','rv']

    for i,var in enumerate(vars):
        for t,file in enumerate(filelist):
            dat = read_data(file,var)
            plot2d(T,X,Z,t,dat,var,folder)
                
    for t,file in enumerate(filelist):
        th = read_data(file,'th')
        u = read_data(file,'u')
        w = read_data(file,'w')
        plot2dquiver(T,X,Z,t,th,u,w,folder)

    #make_gif(folder)

####
# functions
####

# 2d xz profile plots
def plot2d(t,x,z,time,var,name,folder):
    # get x and z
    x = x[0:np.shape(var)[0],0:np.shape(var)[1]]
    z = z[0:np.shape(var)[0],0:np.shape(var)[1]]

    # define plot
    fig = plt.figure(figsize=(8,6))
    fs = 15

    # define dictionaries for plotting
    colors = {'cwc':'Reds', 'rwc':'Blues', 'nc':'Reds', 'nr':'Blues', 'ef':'Purples', 'rv':'Greens'}
    varnames = {'cwc':'cloud water mixing ratio [g/kg]',
        'rwc':'rain water mixing ratio [g/kg]',
        'nc':'cloud droplet concentration [mg$^{-1}$]',
        'nr':'rain droplet concentration [mg$^{-1}$]',
        'ef':'cloud droplet effective radius [$\\mu$m]',
        'rv':'water vapor mixing ratio [g/kg]'}
    #varmins = {'cwc':0, 'rwc':0, 'nc':0, 'nr':0.1, 'ef':0, 'rv':0}
    #varmaxs = {'cwc':2, 'rwc':0.5, 'nc':350, 'nr':10, 'ef':30, 'rv':15}

    varmins = {'cwc':0, 'rwc':0, 'nc':0, 'nr':0.1, 'ef':0, 'rv':0}
    varmaxs = {'cwc':2, 'rwc':0.5, 'nc':200, 'nr':10, 'ef':30, 'rv':20}
    
    # extract from dictionaries
    color = colors[name]
    title = varnames[name]
    vmin = varmins[name]
    vmax = varmaxs[name]

    # set colors
    cmap = mpl.cm.gnuplot2
    color = mpl.colors.LinearSegmentedColormap('gnuplot2', mpl.cm.revcmap(cmap._segmentdata))

    # plot
    if name == 'nr':
        levs = np.geomspace(vmin,vmax,40,endpoint=True)
    else:
        levs = np.linspace(vmin,vmax,40,endpoint=True)
    pc = plt.contourf(x,z,var,levels=levs,cmap=color)

    # set plot attributes
    plt.xlabel('x [m]',fontsize=fs)
    plt.ylabel('z [m]',fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title(title+' , t = '+str(t[time]), fontsize=12)
    
    # make colobar
    if name == 'nr':
        tick = np.geomspace(vmin,vmax,5,endpoint=True)    
    else:
        tick=np.linspace(vmin,vmax,5,endpoint=True)
    plt.colorbar(pc,ticks=tick)

    # save figure
    plt.tight_layout()
    if not os.path.exists(folder+'plots/'+name+'/'):
        os.makedirs(folder+'plots/'+name+'/')
    plt.savefig(folder+'plots/'+name+'/'+str(int(round(t[time])))+'.png',dpi=300)
    plt.close()

# 2dquiver xz profile plot of theta and velocity
def plot2dquiver(t,x,z,time,var,u,w,folder):
    name = 'th_vel'

    # get x and z
    x = x[0:np.shape(var)[0],0:np.shape(var)[1]]
    z = z[0:np.shape(var)[0],0:np.shape(var)[1]]

    # define plot
    fig = plt.figure(figsize=(8,6))
    fs = 15

    # extract from dictionaries
    color = 'Reds'
    title = '$\\theta$, potential temperature (K) and velocity'
    vmin = 285.0
    vmax = 305.0

    # set colors
    cmap = mpl.cm.gnuplot2
    color = mpl.colors.LinearSegmentedColormap('gnuplot2', mpl.cm.revcmap(cmap._segmentdata))

    # contour plot of theta
    levs = np.linspace(vmin,vmax,40,endpoint=True)
    pc = plt.contourf(x,z,var,levels=levs,cmap=color)

    # quiver plot of velocity
    plt.quiver(x,z,u,w)

    # set plot attributes
    plt.xlabel('x [m]',fontsize=fs)
    plt.ylabel('z [m]',fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.title(title+' , t = '+str(t[time]), fontsize=12)

    # colorbar
    plt.colorbar(pc,ticks=np.linspace(vmin,vmax,5,endpoint=True))

    # save figure
    plt.tight_layout()
    if not os.path.exists(folder+'plots/'+name+'/'):
        os.makedirs(folder+'plots/'+name+'/')
    plt.savefig(folder+'plots/'+name+'/'+str(int(round(t[time])))+'.png',dpi=300)
    plt.close()

# make gif
def make_gif(folder):
    for var in ['cwc','rwc','nc','nr','ef','rv']:
        gif_name = 'GIF_'+var
        fps = 2
        file_list = glob.glob(folder+'plots/'+var+'/*.png') # Get all the pngs in the current directory
        list.sort(file_list, key=lambda x: int(x.split('_')[-1].split('.png')[0]))
        clip = mpy.ImageSequenceClip(file_list, fps=fps)
        clip.write_gif('./'+folder+'plots/'+var+'/{}.gif'.format(gif_name), fps=fps)

if __name__ == "__main__":
    main()