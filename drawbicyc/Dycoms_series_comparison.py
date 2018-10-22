#import h5py
from scipy.io import netcdf
import numpy as np
from sys import argv
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from bisect import bisect_left
from matplotlib import rc

# activate latex text rendering
rc('text', usetex=True)

def read_my_array(file_obj):
  file_obj.readline() # discarded line with size of the array
  line = file_obj.readline()
  line = line.split(" ")
  del line[0]
  del line[len(line)-1]
  arr = map(float,line)
  return np.array(arr)



def plot_my_array(axarr, plot_iter, time, val, xlabel=None, ylabel=None, varlabel=None, linestyle='-', dashes=(5,2)):
  #, linestyle=":"):
  x = int(plot_iter / nploty)
  y = plot_iter % nploty
  axarr[x, y].set_xlim([0,6])
  axarr[x, y].plot(time, val, label=varlabel, linestyle=linestyle, linewidth=1, dashes=dashes)
  if xlabel:
    axarr[x, y].set_xlabel(xlabel)
  if ylabel:
    axarr[x, y].set_ylabel(ylabel)


def plot_series(var, empty_plots, plot_iter, nplotx, nploty, show_bin=False):

  # read dycoms results
  dycoms_file = netcdf.netcdf_file("DYCOMS_RF02_results/BLCWG_DYCOMS-II_RF02.scalars.nc", "r")
  time = dycoms_file.variables["time"][:,1,1,:].copy() 
  ntime = dycoms_file.variables["ntime"][:,1,1].copy()
  
  groups = np.arange(14)
  itime = np.arange(0, 21600, 60)
  # time in hours
  itime_h = itime / 3600.
  ivar_arr = np.ndarray(shape=(14))
  # mean val
  mvar_arr = np.ndarray(shape=(len(itime)))
  # extrema
  minvar_arr = np.ndarray(shape=(len(itime)))
  maxvar_arr = np.ndarray(shape=(len(itime)))
  # middle two quartiles
  q1var_arr = np.ndarray(shape=(len(itime)))
  q3var_arr = np.ndarray(shape=(len(itime)))
  
  # positions of bin models
  DHARMA_it = 3 # DHARMA_BO, with coalescence efficiency = 1
  RAMS_it = 7
  var_arr = dycoms_file.variables[var][:,1,1,:].copy()

  # calc entrainment rate
  if var == "zi":
    er = var_arr.copy()
    for g in groups:
      if var_arr[g, 0] > 1e35:
        continue
      er[g, 1:ntime[g]-1] = (var_arr[g, 2:ntime[g]] - var_arr[g, 0:ntime[g]-2]) / (time[g, 2:ntime[g]] - time[g, 0:ntime[g]-2])
      er[g, 0] = (var_arr[g, 1] - var_arr[g, 0]) / (time[g, 1] - time[g, 0])
      er[g, ntime[g]-1] = (var_arr[g, ntime[g]-1] - var_arr[g, ntime[g]-2]) / (time[g, ntime[g]-1] - time[g, ntime[g]-2])
      var_arr[g, 0:ntime[g]] = (er[g,0:ntime[g]] + 3.75e-6 * var_arr[g,0:ntime[g]]) * 100# add LS subsidence and change to cm

  # surf precip - change from W/m2 to mm/d
  rhow = 1e3 # kg/m3
  Lc = 2264.7e3 # J/kg
  if var == "precip":
    var_arr = var_arr / (rhow * Lc) * 1e3 * 24 * 3600
  
  # interpolate to same time positions
  mean_iter = 0
  for it in itime:
    for g in groups:
#      if var_arr[g, 0] > 1e35:  #netcdf fill values are read as ca. 9e36
#        continue
      # find index with time less than it
      i = bisect_left(time[g, 0:ntime[g]], it)
      # same time
      if time[g, i] == it:
        ivar_arr[g] = var_arr[g, i]
      # time[g, i] > it
      else:
        if i == 0:
          ivar_arr[g] = var_arr[g,i]
        else:
          prev_time = time[g, i-1]
          prev_var_arr = var_arr[g, i-1]
          ivar_arr[g] = prev_var_arr + (it - prev_time) / (time[g, i] - prev_time) * (var_arr[g, i] - prev_var_arr)
    mvar_arr[mean_iter] = ivar_arr[ivar_arr < 1e35].mean() # < 1e35 to avoid the netcdf fill values from models that didn't calculate this vat
    minvar_arr[mean_iter] = ivar_arr[ivar_arr < 1e35].min()
    maxvar_arr[mean_iter] = ivar_arr[ivar_arr < 1e35].max()
    q1var_arr[mean_iter] = np.percentile(ivar_arr[ivar_arr < 1e35], 25)
    q3var_arr[mean_iter] = np.percentile(ivar_arr[ivar_arr < 1e35], 75)
    mean_iter+=1

  x = int(plot_iter / nploty)
  y = plot_iter % nploty
  
  #  for g in groups:
  #    print var
  #    print var_arr[g, 0:ntime[g]]
  #    if var_arr[g,0] < 1e35: # netcdf fill values are erad as ca. 9e36
  #      axarr[x, y].plot(time[g,0:ntime[g]] / 3600., var_arr[g,0:ntime[g]])
    
  axarr[x, y].fill_between(itime_h, minvar_arr, maxvar_arr, color='0.9')
  axarr[x, y].fill_between(itime_h, q1var_arr, q3var_arr, color='0.7')
  axarr[x, y].plot(itime_h, mvar_arr, color='black')
  # plot precip and NC of bin models
  if show_bin:
    if False:#True:# var == "precip" or var == "ndrop_cld":
      DHARMA_time = time[DHARMA_it,0:ntime[DHARMA_it]].copy() / 3600.
      DHARMA_precip = var_arr[DHARMA_it,0:ntime[DHARMA_it]].copy()
   #   print DHARMA_time, DHARMA_precip
      axarr[x, y].plot(DHARMA_time[:], DHARMA_precip[:], color='red', linewidth=1)
      RAMS_time = time[RAMS_it,0:ntime[RAMS_it]].copy() / 3600.
      RAMS_precip = var_arr[RAMS_it,0:ntime[RAMS_it]].copy()
    #  print RAMS_time, RAMS_precip
      axarr[x, y].plot(RAMS_time[:], RAMS_precip[:], color='green', linewidth=1)

  
  dycoms_file.close()
  
  
  #read my results
  series_files_names = []
  series_labels = []
  file_no = np.arange(1, len(sys.argv)-1 , 2)
  for no in file_no:
    series_files_names.append(argv[no])
    series_labels.append(argv[no+1])
  
  label_counter=0
  for file_name in series_files_names:
    
    series_file = open(file_name, "r")
    my_times = read_my_array(series_file)
    my_max_w_var = read_my_array(series_file)
    my_cfrac = read_my_array(series_file)
    my_lwp = read_my_array(series_file)
    my_er = read_my_array(series_file)
    my_sp = read_my_array(series_file)
    read_my_array(series_file) # discard accumulated precip
    my_act_cond = read_my_array(series_file)
    my_zb = read_my_array(series_file)
    
    # rescale time to hours
    my_times = my_times / 3600.
    
    series_file.close()
    
    linestyles = ['--', '-.', ':']
    dashList = [(3,1),(1,1),(4,1,1,1),(4,2)] 
    if var == "lwp":
      plot_my_array(axarr, plot_iter, my_times, my_lwp, ylabel='LWP [g m$^{-2}$]', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
    if var == "zi":
      plot_my_array(axarr, plot_iter, my_times, my_er, ylabel='Entrainment rate [cm s$^{-1}$]', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
    if var == "w2_max":
      plot_my_array(axarr, plot_iter, my_times, my_max_w_var, ylabel='Max. $w$ variance [m$^{2}$ s$^{-2}$]', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
    if var == "precip":
      plot_my_array(axarr, plot_iter, my_times, my_sp, xlabel='Time [h]', ylabel='Surface precip. [mm / day]', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
    if var == "ndrop_cld":
      plot_my_array(axarr, plot_iter, my_times, my_act_cond, xlabel='Time [h]', ylabel='$N_c$ [cm$^{-3}$]', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
    if var == "zb":
      plot_my_array(axarr, plot_iter, my_times, my_zb, xlabel='Time [h]', ylabel='Cloud base height [m]', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
  # plot_my_array(axarr, plot_iter, my_times, my_cfrac, ylabel='Cloud fraction', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
    label_counter+=1

  plot_iter = plot_iter + 1
  return plot_iter



# init the plot
nplotx = 2
nploty= 3
fig, axarr = plt.subplots(nplotx,nploty)


dycoms_vars = ["lwp", "zi", "w2_max", "precip", "ndrop_cld", "zb"]# "cfrac"]

if len(dycoms_vars) % nploty == 0:
  nemptyplots = 0
else:
  nemptyplots = nploty - len(dycoms_vars) % nploty
emptyplots = np.arange(nploty - nemptyplots, nploty)

plot_iter = 0

for var in dycoms_vars:
  plot_iter = plot_series(var, emptyplots, plot_iter,nplotx, nploty)
  print plot_iter


# show legends on each subplot
#for x in np.arange(nplotx):
#  for y in np.arange(nploty):
#    axarr[x,y].legend()

# legend font size
plt.rcParams.update({'font.size': 8})

# hide axes on empty plots
for empty in emptyplots:
  axarr[nplotx-1, empty].axis('off')

# hide hrzntl tic labels
x_empty_label = np.arange(0, nplotx-1)
y_empty_label = np.arange(nploty)
for x in x_empty_label:
  for y in y_empty_label:
    axarr[x,y].set_xticklabels([])

#axes = plt.gca()
#axes.tick_params(direction='in')
x_arr = np.arange(nplotx)
y_arr = np.arange(nploty)
for x in x_arr:
  for y in y_arr:
    #tics inside
    axarr[x,y].tick_params(direction='in', which='both', top=1, right=1)
    #minor tics
    axarr[x,y].xaxis.set_minor_locator(AutoMinorLocator())
    axarr[x,y].yaxis.set_minor_locator(AutoMinorLocator())
    #labels and tics font size
    for item in ([axarr[x,y].xaxis.label, axarr[x,y].yaxis.label] + axarr[x,y].get_xticklabels() + axarr[x,y].get_yticklabels()):
      item.set_fontsize(8)

#single legend for the whole figure
handles, labels = axarr[0,0].get_legend_handles_labels()

lgd = fig.legend(handles, labels, handlelength=4, loc='lower center', bbox_to_anchor=(0.45,0))

#figure size
fig.set_size_inches(7.874, 5. + (len(labels) - 2) * 0.2)# 5.214)#20.75,13.74)
#distances between subplots and from bottom of the plot
fig.subplots_adjust(bottom=0.18 + (len(labels) - 2) * 0.03, hspace=0.1, wspace=0.4)

#fig.tight_layout(pad=0, w_pad=0, h_pad=0)

#plt.show()
fig.savefig(argv[len(sys.argv)-1], bbox_inches='tight', dpi=300)#, bbox_extra_artists=(lgd,))
