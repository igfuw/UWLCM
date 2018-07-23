#import h5py
from scipy.io import netcdf
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from bisect import bisect_left


def read_my_array(file_obj):
  series_file.readline() # discarded line with size of the array
  line = series_file.readline()
  line = line.split(" ")
  del line[0]
  del line[len(line)-1]
  arr = map(float,line)
  return np.array(arr)

def plot_my_array(axarr, plot_iter, time, val, xlabel=None, ylabel=None ):
  x = int(plot_iter / 2)
  y = plot_iter % 2
  axarr[x, y].plot(time, val, color='red')
  if xlabel:
    axarr[x, y].set_xlabel(xlabel)
  if ylabel:
    axarr[x, y].set_ylabel(ylabel)
  plot_iter = plot_iter+1
  return plot_iter



# init the plot
fig, axarr = plt.subplots(3,2)


# read dycoms results

dycoms_file = netcdf.netcdf_file("DYCOMS_RF02_results/BLCWG_DYCOMS-II_RF02.scalars.nc", "r")
time = dycoms_file.variables["time"][:,1,1,:].copy() 
ntime = dycoms_file.variables["ntime"][:,1,1].copy()

dycoms_vars = ["lwp", "w2_max", "precip", "ndrop_cld", "cfrac", "zi"]

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

plot_iter = 0
for var in dycoms_vars:
  var_arr = dycoms_file.variables[var][:,1,1,:].copy()

  # calc entrainment rate
  if var == "zi":
    er = var_arr.copy()
    for g in groups:
      er[g, 1:ntime[g]-1] = (var_arr[g, 2:ntime[g]] - var_arr[g, 0:ntime[g]-2]) / (time[g, 2:ntime[g]] - time[g, 0:ntime[g]-2])
      er[g, 0] = (var_arr[g, 1] - var_arr[g, 0]) / (time[g, 1] - time[g, 0])
      er[g, ntime[g]-1] = (var_arr[g, ntime[g]-1] - var_arr[g, ntime[g]-2]) / (time[g, ntime[g]-1] - time[g, ntime[g]-2])
    var_arr = (er + 3.75e-6 * var_arr) * 100# add LS subsidence and change to cm

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
    mvar_arr[mean_iter] = ivar_arr[ivar_arr < 1e35].mean()
    minvar_arr[mean_iter] = ivar_arr[ivar_arr < 1e35].min()
    maxvar_arr[mean_iter] = ivar_arr[ivar_arr < 1e35].max()
    q1var_arr[mean_iter] = np.percentile(ivar_arr[ivar_arr < 1e35], 25)
    q3var_arr[mean_iter] = np.percentile(ivar_arr[ivar_arr < 1e35], 75)
    mean_iter+=1

  x = int(plot_iter / 2)
  y = plot_iter % 2

#  for g in groups:
#    print var
#    print var_arr[g, 0:ntime[g]]
#    if var_arr[g,0] < 1e35: # netcdf fill values are erad as ca. 9e36
#      axarr[x, y].plot(time[g,0:ntime[g]] / 3600., var_arr[g,0:ntime[g]])
  
  axarr[x, y].fill_between(itime_h, minvar_arr, maxvar_arr, color='0.9')
  axarr[x, y].fill_between(itime_h, q1var_arr, q3var_arr, color='0.7')
  axarr[x, y].plot(itime_h, mvar_arr, color='black')

  plot_iter += 1

dycoms_file.close()


#read my results

series_file = open(argv[1], "r")
my_times = read_my_array(series_file)
my_max_w_var = read_my_array(series_file)
my_cfrac = read_my_array(series_file)
my_lwp = read_my_array(series_file)
my_er = read_my_array(series_file)
my_sp = read_my_array(series_file)
read_my_array(series_file) # discard accumulated precip
my_act_cond = read_my_array(series_file)

# rescale time to hours
my_times = my_times / 3600.

series_file.close()


#dycoms_vars = ["lwp", "w2_max", "precip", "ndrop_cld", "cfrac", "zi"]
plot_iter=0
plot_iter = plot_my_array(axarr, plot_iter, my_times, my_lwp, ylabel='LWP (g / m$^{2}$)')
plot_iter = plot_my_array(axarr, plot_iter, my_times, my_max_w_var, ylabel='max w variance (m$^{2}$ / s$^2$)')
plot_iter = plot_my_array(axarr, plot_iter, my_times, my_sp, ylabel='surface precip. (mm / day)')
plot_iter = plot_my_array(axarr, plot_iter, my_times, my_act_cond, ylabel='N$_c$ (cm$^{-3}$)')
plot_iter = plot_my_array(axarr, plot_iter, my_times, my_cfrac, xlabel='time (h)', ylabel='cloud frac.')
plot_iter = plot_my_array(axarr, plot_iter, my_times, my_er, xlabel='time (h)', ylabel='entrainment rate (cm / s)')


#axarr[0, 0].set_xlabel('time(h)')


plt.show()

#nx, ny, nz = w3d.shape

#for t in range(time_start, time_end+1, outfreq):
#  filename = directory + "/timestep" + str(t).zfill(10) + ".h5"
#  print filename
#
#  for vel in velocities:
#
#    w3d = h5py.File(filename, "r")[vel][:,:,:]
#    
#    
#    w2d = w3d[:, :, poziom]
#    
#    wkx = 1.0 / np.sqrt(nx - 1) * np.fft.rfft(w2d, axis = 0)
#    wky = 1.0 / np.sqrt(ny - 1) * np.fft.rfft(w2d, axis = 1)
#    
#    Ex = 0.5 * (np.abs(wkx) ** 2)
#    Ex = np.mean(Ex, axis = 1)
#    Ey = 0.5 * (np.abs(wky) ** 2)
#    Ey = np.mean(Ey, axis = 0)
#    
#    Exy = 0.5 * (Ex + Ey)
#    Exy_avg[vel] += Exy
#
#    K = np.fft.rfftfreq(nx - 1)
##    plt.loglog(K, Exy)
#  
#  if (t == time_start):
#    plt.loglog(K, 1e-4*K**(-5./3.))
#
#for vel in velocities:
#  Exy_avg[vel] /= (time_end - time_start) / outfreq + 1
#  plt.loglog(K, Exy_avg[vel], linewidth=2)
# 
#plt.show()
