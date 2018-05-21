#import h5py
from scipy.io import netcdf
import numpy as np
from sys import argv
import sys
import matplotlib.pyplot as plt
from bisect import bisect_left


def read_my_array(file_obj):
  file_obj.readline() # discarded line with size of the array
  line = file_obj.readline()
  line = line.split(" ")
  del line[0]
  del line[len(line)-1]
  arr = map(float,line)
  return np.array(arr)



dycoms_vars = ["thetal", "qt", "ql", "w_var", "w_skw", "precip", "ss", "cfrac", "ndrop_cld"]

# init the plot
nplotx = 2 #int(len(dycoms_vars)/5 + 1)
nploty = int(float(len(dycoms_vars))/float(nplotx) + 0.5)
nemptyplots = nploty - len(dycoms_vars) % nploty
emptyplots = np.arange(nploty - nemptyplots, nploty)
print nplotx, nploty, nemptyplots
print emptyplots
fig, axarr = plt.subplots(nplotx, nploty )

def plot_my_array(axarr, plot_iter, time, val, xlabel=None, ylabel=None, varlabel=None ):
  x = int(plot_iter / nploty)
  y = plot_iter % nploty
  if varlabel != None:
    axarr[x, y].plot(time, val, label=varlabel)
  else:
    axarr[x, y].plot(time, val)
  if xlabel:
    axarr[x, y].set_xlabel(xlabel)
  if ylabel:
    axarr[x, y].set_ylabel(ylabel)
  plot_iter = plot_iter+1
  return plot_iter


# read dycoms results

dycoms_file = netcdf.netcdf_file("DYCOMS_RF02_results/BLCWG_DYCOMS-II_RF02.profiles.nc", "r")
dycoms_series_file = netcdf.netcdf_file("DYCOMS_RF02_results/BLCWG_DYCOMS-II_RF02.scalars.nc", "r")
time = dycoms_file.variables["time"][:].copy() 
zt = dycoms_file.variables["zt"][:,1,1,:].copy() 
nzt = dycoms_file.variables["nzt"][:,1,1].copy()

series_time = dycoms_series_file.variables["time"][:, 1, 1, :].copy() 
series_zi = dycoms_series_file.variables["zi"][:, 1, 1, :].copy() 
series_ntime = dycoms_series_file.variables["ntime"][:, 1, 1].copy() 

ntime = 13

# at each time, zt needs to be rescaled by inversion height, this rescaled value will be stored here
rzt = np.zeros((301)) # group idx/ height idx


groups = np.arange(14)
ihght = np.arange(0, 1.2, 0.01) # height levels scaled by inversionb height to which we will inteprolate results

ivar_arr = np.ndarray(shape=(14, len(ihght))) # to store interpolated average over time, group idx / height idx
#ivar_arr = np.ndarray(shape=(14, 13)) # group index / time index

# mean val
mvar_arr = np.ndarray(shape=(len(ihght)))
# extrema
minvar_arr = np.ndarray(shape=(len(ihght)))
maxvar_arr = np.ndarray(shape=(len(ihght)))
# middle two quartiles
q1var_arr = np.ndarray(shape=(len(ihght)))
q3var_arr = np.ndarray(shape=(len(ihght)))

plot_iter = 0
for var in dycoms_vars:
  var_arr = dycoms_file.variables[var][:,1,1,:,:].copy() 

  for g in groups:
    ivar_arr[g,:] = 0.
    time_index = 0
    for t in time:
      # we start averaging after 2h
      if t <= 7200.:
        time_index += 1
        continue
      # -- rescale height by inversion height at fiven time --
      # initial iinversion height is 795m
      if t==0:
        rzt[:] = zt[g, :] / 795
      else:
        # find nearest zi from time series, TODO: we should interpolate here
        series_time_index = bisect_left(series_time[g, 0:series_ntime[g]], t)
        #series_time_index = np.where(series_time[g,0:series_ntime[g]]==t) 
        rzt[:] = zt[g, :] / series_zi[g, series_time_index]

      # interpolate to same height positions and add to the mean over time fir this group
      i_hght_index = 0
      for it in ihght:
        # find index with height less than it
        rzt_hght_idx = bisect_left(rzt[:], it)
        # same hght
        if rzt[rzt_hght_idx] == it:
          ivar_arr[g, i_hght_index] += var_arr[g, time_index, rzt_hght_idx]
        # time[g, i] > it
        else:
          if rzt_hght_idx == 0:
            ivar_arr[g, i_hght_index] += var_arr[g, time_index, rzt_hght_idx]
          else:
            prev_hght = rzt[rzt_hght_idx-1]
            prev_var_arr = var_arr[g, time_index, rzt_hght_idx-1]
            ivar_arr[g, i_hght_index] +=  prev_var_arr + (it - prev_hght) / (rzt[rzt_hght_idx] - prev_hght) * (var_arr[g, time_index, rzt_hght_idx]- prev_var_arr)
        i_hght_index += 1
      time_index += 1
    ivar_arr[g,:] /= (13-5) 
#    axarr[0, 0].plot(ivar_arr[g,:], ihght)

  # calc statistics from groups

 # mean_iter = 0
  #mvar_arr[mean_iter] = ivar_arr[ivar_arr < 1e35].mean() # < 1e35 to avoid the netcdf fill values from models that didn't calculate this vat
  for zi in np.arange(len(ihght)):
    ivar_arr_1d = ivar_arr[:,zi]
    mvar_arr[zi] = ivar_arr_1d[ivar_arr_1d < 1e35].mean() # < 1e35 to avoid the netcdf fill values from models that didn't calculate this vat
    minvar_arr[zi] = ivar_arr_1d[ivar_arr_1d < 1e35].min()
    maxvar_arr[zi] = ivar_arr_1d[ivar_arr_1d < 1e35].max()
    q1var_arr[zi] = np.percentile(ivar_arr_1d[ivar_arr_1d < 1e35], 25)
    q3var_arr[zi] = np.percentile(ivar_arr_1d[ivar_arr_1d < 1e35], 75)
  #mvar_arr[zi] = np.mean(ivar_arr[:,:], axis=0)
#  minvar_arr[:] = np.min(ivar_arr[:,:], axis=0)
#  maxvar_arr[:] = np.max(ivar_arr[:,:], axis=0)
#  q1var_arr[:] = np.percentile(ivar_arr[:,:], 25, axis=0)
#  q3var_arr[:] = np.percentile(ivar_arr[:,:], 75, axis=0)
 # mean_iter+=1

  x = int(plot_iter / nploty)
  y = plot_iter % nploty

#  for g in groups:
#    print var
#    print var_arr[g, 0:ntime[g]]
#    if var_arr[g,0] < 1e35: # netcdf fill values are erad as ca. 9e36
#      axarr[x, y].plot(time[g,0:ntime[g]] / 3600., var_arr[g,0:ntime[g]])
  
  axarr[x, y].fill_betweenx(ihght, minvar_arr, maxvar_arr, color='0.9')
  axarr[x, y].fill_betweenx(ihght, q1var_arr, q3var_arr, color='0.7')
  axarr[x, y].plot(mvar_arr, ihght, color='black')
#  axes = axarr[x, y].gca()
  axarr[x, y].set_ylim([0,1.2])
  if var == "ss":
    #axarr[x, y].set_xlim(xmin=-1)
    axarr[x, y].set_xlim([-2,1])
#
#  for g in groups:
#    axarr[x, y].plot(mvar_arr[g])#[mvar_arr[g]<1e35])

  plot_iter += 1

dycoms_file.close()


#read my results
profiles_files_names = []
file_no = np.arange(len(sys.argv)-1)
for no in file_no:
  profiles_files_names.append(argv[no+1])

label_counter = 0
for file_name in profiles_files_names:
  
  # dycoms_vars = ["thetal", "qt", "ql", "w_var", "w_skw", "precip", "ss", "cfrac", "ndrop_cld", "qr"]
  try:
    profiles_file = open(file_name, "r")
    my_pos = read_my_array(profiles_file)
    my_rtot = read_my_array(profiles_file)
    my_rliq = read_my_array(profiles_file)
    my_thl = read_my_array(profiles_file)
    my_wvar = read_my_array(profiles_file)
    my_w3rd = read_my_array(profiles_file)
    my_prflux = read_my_array(profiles_file)
    my_clfrac = read_my_array(profiles_file)
    my_nc = read_my_array(profiles_file)
    my_ss = read_my_array(profiles_file)

    print 'mean nc in cloud cells: ' , np.mean(my_nc[my_nc>20])
  
    profiles_file.close()
    
  
    plot_iter=0
    plot_iter = plot_my_array(axarr, plot_iter, my_thl, my_pos, xlabel='$\theta_l$[K]', varlabel=label_counter)
    plot_iter = plot_my_array(axarr, plot_iter, my_rtot, my_pos, xlabel='q$_{tot}$[g/kg]', varlabel=label_counter)
    plot_iter = plot_my_array(axarr, plot_iter, my_rliq, my_pos, xlabel='q$_{l}$[g/kg]', varlabel=label_counter)
    plot_iter = plot_my_array(axarr, plot_iter, my_wvar, my_pos, xlabel='wvar', varlabel=label_counter)
    plot_iter = plot_my_array(axarr, plot_iter, my_w3rd, my_pos, xlabel='w3rd', varlabel=label_counter)
    plot_iter = plot_my_array(axarr, plot_iter, my_prflux, my_pos, xlabel='prflux', varlabel=label_counter)
    plot_iter = plot_my_array(axarr, plot_iter, my_ss, my_pos, xlabel='S', varlabel=label_counter)
    plot_iter = plot_my_array(axarr, plot_iter, my_clfrac, my_pos, xlabel='cl frac', varlabel=label_counter)
    plot_iter = plot_my_array(axarr, plot_iter, my_nc, my_pos, xlabel='nc', varlabel=label_counter)
  except:
    print 'error opening file: ', file_name
    my_pos = 0
    my_rtot = 0
    my_rliq = 0
    my_thl = 0
    my_wvar = 0
    my_w3rd = 0
    my_prflux = 0
    my_clfrac = 0
    my_nc = 0
    my_ss = 0
  label_counter = label_counter+1

# hide axes on empty plots
for empty in emptyplots:
  axarr[nplotx-1, empty].axis('off')
# show legends
for x in np.arange(nplotx):
  for y in np.arange(nploty):
    axarr[x,y].legend()
plt.show()

