from scipy.io import netcdf
import numpy as np
from sys import argv
from bisect import bisect_left
import os

# labels used in the Dycoms intercomparison reslts
# empty label means that there is no result available
dycoms_labels = {
  "thl" : "thetal",
  "00rtot" : "qt",
  "rliq" : "ql",
  "prflux" : "precip",
  "cl_nc" : "ndrop_cld",
  "clfrac" : "cfrac",
  "wvar" : "w_var",
  "w3rd" : "w_skw",
  "sat_RH" : "ss",
  "rad_flx" : "rad_flx",
  "lwp" : "lwp",
  "er" : "zi",
  "wvarmax" : "w2_max",
  "surf_precip" : "precip",
  "cloud_base" : "zb",
  "cl_gccn_conc" : "",
  "gccn_rw_cl" : "",
  "non_gccn_rw_cl" : ""
}

def plot_reference_profiles(var_list, plot_iter, nplotx, nploty, axarr, show_bin=False):
  # files with Dycoms intercomparison results
  dir_path = os.path.dirname(os.path.realpath(__file__))
  dycoms_file = netcdf.netcdf_file(dir_path+"/DYCOMS_RF02_results/BLCWG_DYCOMS-II_RF02.profiles.nc", "r")
  dycoms_series_file = netcdf.netcdf_file(dir_path+"/DYCOMS_RF02_results/BLCWG_DYCOMS-II_RF02.scalars.nc", "r")
  time = dycoms_file.variables["time"][:].copy() 
  zt = dycoms_file.variables["zt"][:,1,1,:].copy() 
  nzt = dycoms_file.variables["nzt"][:,1,1].copy()
  
  series_time = dycoms_series_file.variables["time"][:, 1, 1, :].copy() 
  series_zi = dycoms_series_file.variables["zi"][:, 1, 1, :].copy() 
  series_ntime = dycoms_series_file.variables["ntime"][:, 1, 1].copy() 
  
  # positions of bin models
  DHARMA_it = 3 # DHARMA_BO, with coalescence efficiency = 1
  RAMS_it = 7
  
  ntime = 13
  
  for var in var_list:
    x = int(plot_iter / nploty)
    y = plot_iter % nploty

    # read dycoms results
    if dycoms_labels[var] != '': # empty dycoms_label indicates that this plot is not available from the dycoms intercomparison
      # at each time, zt needs to be rescaled by inversion height, this rescaled value will be stored here
      rzt = np.zeros((301)) # group idx/ height idx
      
      groups = np.arange(14)
      ihght = np.arange(0, 1.6, 0.01) # height levels scaled by inversionb height to which we will inteprolate results
      
      ivar_arr = np.ndarray(shape=(14, len(ihght))) # to store interpolated average over time, group idx / height idx
      
      # mean val
      mvar_arr = np.ndarray(shape=(len(ihght)))
      # extrema
      minvar_arr = np.ndarray(shape=(len(ihght)))
      maxvar_arr = np.ndarray(shape=(len(ihght)))
      # middle two quartiles
      q1var_arr = np.ndarray(shape=(len(ihght)))
      q3var_arr = np.ndarray(shape=(len(ihght)))
      
      
      if var=="cl_nc_zoom":
        var="cl_nc"
      var_arr = dycoms_file.variables[dycoms_labels[var]][:,1,1,:,:].copy() 
      
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

        # plot precip and NC of bin models
        if show_bin and g == DHARMA_it:
          DHARMA_prof = ivar_arr[g,:]
          DHARMA_pos = ihght[:]
          axarr[x, y].plot(DHARMA_prof[DHARMA_prof < 1e35], DHARMA_pos[DHARMA_pos < 1e35], linewidth=1, label="DHARMA", color='red')
        if show_bin and g == RAMS_it:
          RAMS_prof = ivar_arr[g,:]
          RAMS_pos = ihght[:]
          axarr[x, y].plot(RAMS_prof[RAMS_prof < 1e35], RAMS_pos[RAMS_pos < 1e35], linewidth=1, label="RAMS", color='green')
      
      # calc statistics from groups
      
      for zi in np.arange(len(ihght)):
        ivar_arr_1d = ivar_arr[:,zi]
        mvar_arr[zi] = ivar_arr_1d[ivar_arr_1d < 1e35].mean() # < 1e35 to avoid the netcdf fill values from models that didn't calculate this vat
        minvar_arr[zi] = ivar_arr_1d[ivar_arr_1d < 1e35].min()
        maxvar_arr[zi] = ivar_arr_1d[ivar_arr_1d < 1e35].max()
        q1var_arr[zi] = np.percentile(ivar_arr_1d[ivar_arr_1d < 1e35], 25)
        q3var_arr[zi] = np.percentile(ivar_arr_1d[ivar_arr_1d < 1e35], 75)
      
      axarr[x, y].fill_betweenx(ihght, minvar_arr, maxvar_arr, color='0.9')
      axarr[x, y].fill_betweenx(ihght, q1var_arr, q3var_arr, color='0.7')
      axarr[x, y].plot(mvar_arr, ihght, color='black')
    plot_iter+=1
  dycoms_file.close()

def plot_reference_series(var_list, plot_iter, nplotx, nploty, axarr, show_bin=False):

  # files with Dycoms intercomparison results
  dir_path = os.path.dirname(os.path.realpath(__file__))
  dycoms_file = netcdf.netcdf_file(dir_path+"/DYCOMS_RF02_results/BLCWG_DYCOMS-II_RF02.scalars.nc", "r")
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

  for var in var_list:
    x = int(plot_iter / nploty)
    y = plot_iter % nploty


    # read dycoms results
    if dycoms_labels[var] != '': # empty dycoms_label indicates that this plot is not available from the dycoms intercomparison
      var_arr = dycoms_file.variables[dycoms_labels[var]][:,1,1,:].copy()
  
      # calc entrainment rate
      if var == "er":
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
      if var == "surf_precip":
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
        DHARMA_time = time[DHARMA_it,0:ntime[DHARMA_it]].copy() / 3600.
        DHARMA_precip = var_arr[DHARMA_it,0:ntime[DHARMA_it]].copy()
      #  print DHARMA_time, DHARMA_precip
        axarr[x, y].plot(DHARMA_time[:], DHARMA_precip[:], color='red', linewidth=1, label="DHARMA")
        RAMS_time = time[RAMS_it,0:ntime[RAMS_it]].copy() / 3600.
        RAMS_precip = var_arr[RAMS_it,0:ntime[RAMS_it]].copy()
       # print RAMS_time, RAMS_precip
        axarr[x, y].plot(RAMS_time[:], RAMS_precip[:], color='green', linewidth=1, label="RAMS")
    plot_iter = plot_iter + 1
  dycoms_file.close()
