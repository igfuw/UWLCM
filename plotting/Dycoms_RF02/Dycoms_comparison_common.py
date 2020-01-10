from scipy.io import netcdf
import numpy as np
from sys import argv
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from bisect import bisect_left
from matplotlib import rc
import os

import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../latex_labels/")
from latex_labels import *

# labels used in the Dycoms intercomparison reslts
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

labeldict = {
 0 : "(a)",
 1 : "(b)",
 2 : "(c)",
 3 : "(d)",
 4 : "(e)",
 5 : "(f)",
 6 : "(g)",
 7 : "(h)",
 8 : "(i)",
 9 : "(j)"
}

yscaledict = {
  "thl" : "linear",
  "00rtot" : "linear",
  "rliq" : "linear",
  "prflux" : "linear",
  "cl_nc" : "linear",
  "clfrac" : "linear",
  "wvar" : "linear",
  "w3rd" : "linear",
  "sat_RH" : "linear",
  "rad_flx" : "linear",
  "lwp" : "linear",
  "er" : "linear",
  "wvarmax" : "linear",
  "surf_precip" : "linear",
  "acc_precip" : "linear",
  "cloud_base" : "linear",
  "gccn_rw_cl" : "linear",
  "non_gccn_rw_cl" : "linear",
  "base_prflux_vs_clhght" : "log",
  "cl_gccn_conc" : "log"
}

xlimdict_profs = {
  "thl" : None,
  "00rtot" : None,
  "rliq" : None,
  "prflux" : (0,20),
  "cl_nc" : (0,90),
  "clfrac" : None,
  "wvar" : None,
  "w3rd" : None,
  "sat_RH" : None,
  "rad_flx" : None,
  "gccn_rw_cl" : (0,90),
  "non_gccn_rw_cl" : (0,12),
  "base_prflux_vs_clhght" : (1,10000)
}

ylimdict_profs = {
  "thl" : (0,3000),
  "00rtot" : (0,3000),
  "rliq" : (0,3000),
  "prflux" : (0,3000),
  "cl_nc" : (0,3000),
  "clfrac" : (0,3000),
  "wvar" : (0,3000),
  "w3rd" : (0,3000),
  "sat_RH" : (0,3000),
  "rad_flx" : (0,3000),
  "gccn_rw_cl" : (0,3000),
  "non_gccn_rw_cl" : (0,3000),
  "base_prflux_vs_clhght" : (0,2500)
}

xlimdict_series = {
  "clfrac" : (1,5),
  "lwp" : (1,5),
  "er" : (1,5),
  "wvarmax" : (1,5),
  "surf_precip" : (1,5),
  "acc_precip" : (1,5),
  "cloud_base" : (1,5),
  "cl_gccn_conc" : (1,5)
}

ylimdict_series = {
  "clfrac" : None,
  "lwp" : None,
  "er" : None,
  "wvarmax" : None,
  "surf_precip" : None,
  "acc_precip" : (0,0.07),
  "cloud_base" : None,
  "cl_gccn_conc" : None
}

def read_my_array(file_obj):
  arr_name = file_obj.readline()
  file_obj.readline() # discarded line with size of the array
  line = file_obj.readline()
  line = line.split(" ")
  del line[0]
  del line[len(line)-1]
  arr = map(float,line)
  return np.array(arr), arr_name

def read_my_var(file_obj, var_name):
  while True:
    arr, name = read_my_array(file_obj)
    if(str(name).strip() == str(var_name).strip()):
      break
  return arr

def plot_my_array(axarr, plot_iter, time, val, nploty, xlabel=None, ylabel=None, varlabel=None , linestyle='--', dashes=(5,2), xlim=None, ylim=None, xscale="linear"):
  x = int(plot_iter / nploty)
  y = plot_iter % nploty
  if varlabel != None:
    axarr[x, y].plot(time, val, label=varlabel, linestyle=linestyle, linewidth=1, dashes=dashes)
  else:
    axarr[x, y].plot(time, val, linestyle=linestyle, linewidth=1, dashes=dashes)
  if xlabel:
    axarr[x, y].set_xlabel(xlabel)
  if ylabel:
    axarr[x, y].set_ylabel(ylabel)
  if xlim:
    axarr[x, y].set_xlim(xlim)
  if ylim:
    axarr[x, y].set_ylim(ylim)
  axarr[x, y].set_xscale(xscale)

def plot_profiles(var_list, plot_iter, nplotx, nploty, axarr, show_bin=False, suffix='', reference=True, ylabel=''):
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

  # create a list of files with UWLCM results
  profiles_files_names = []
  profiles_labels = []
  file_no = np.arange(1, len(sys.argv)-2 , 2)
  print file_no
  for no in file_no:
    profiles_files_names.append(argv[no]+suffix)
    profiles_labels.append(argv[no+1])
  
  for var in var_list:
    x = int(plot_iter / nploty)
    y = plot_iter % nploty

    # set axis ranges
    #if var == "thetal":
    #  axarr[x, y].set_xlim([288.2,289.2])
    #if var == "qt":
    #  axarr[x, y].set_xlim([9,10.4])
    if var == "rliq":
      axarr[x, y].set_xlim([0,.7])
    if var == "wvar":
      axarr[x, y].set_xlim([-0.1,.8])
    if var == "w3rd":
      axarr[x, y].set_xlim([-0.15,.15])
    if var == "sat_RH":
      axarr[x, y].set_xlim([-5,1])
    if var=="cl_nc_zoom":
      axarr[x, y].set_xlim([50,80])
    if var=="cl_nc":
      axarr[x, y].set_xlim([0,120])
    if var=="prflux":
      axarr[x, y].set_xlim([0,70])
    if var=="non_gccn_rw_cl":
      axarr[x, y].set_xlim([0,7])
    if var=="gccn_rw_cl":
      axarr[x, y].set_xlim([0,40])

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
      
    #  if reference:
    #    axarr[x, y].fill_betweenx(ihght, minvar_arr, maxvar_arr, color='0.9')
    #    axarr[x, y].fill_betweenx(ihght, q1var_arr, q3var_arr, color='0.7')
    #    axarr[x, y].plot(mvar_arr, ihght, color='black')
    
    
    axarr[x, y].set_ylim([0,1.2])

    #read my results
    label_counter = 0
    for file_name in profiles_files_names:
    #  try:
      print file_name
      profiles_file = open(file_name, "r")
      my_pos = read_my_var(profiles_file, "position")
      my_res = read_my_var(profiles_file, var)

   #     print 'mean nc in cloud cells: ' , np.mean(my_nc[my_nc>20])

      profiles_file.close()

      linestyles = ['--', '-.', ':']
      dashList = [(3,1),(1,1),(4,1,1,1),(4,2)]

      # get rid of false signal due to cloudiness mask used based on ql instead of nc
      if var == "gccn_rw_cl" or var == "non_gccn_rw_cl":
        my_res[0:3] = 0  

      if(var == "base_prflux_vs_clhght"):
        plot_my_array(axarr, plot_iter, my_res, my_pos, nploty, xlabel=var_labels[var], ylabel="cloudy column height [m]", varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)], xlim = (0,600))
      else:
        plot_my_array(axarr, plot_iter, my_res, my_pos, nploty, xlabel=var_labels[var], ylabel=ylabel, varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
   #   except:
   #     print 'error opening file: ', file_name
   #     my_pos = 0
   #     my_res = 0
      label_counter = label_counter+1

    plot_iter += 1

  
#  label_counter = 0
#  for file_name in profiles_files_names:
#    
#    try:
#      profiles_file = open(file_name, "r")
#      my_pos = read_my_array(profiles_file)
#      my_rtot = read_my_array(profiles_file)
#      my_rliq = read_my_array(profiles_file)
#      my_thl = read_my_array(profiles_file)
#      my_wvar = read_my_array(profiles_file)
#      my_w3rd = read_my_array(profiles_file)
#      my_prflux = read_my_array(profiles_file)
#      my_clfrac = read_my_array(profiles_file)
#      my_nc = read_my_array(profiles_file)
#      my_ss = read_my_array(profiles_file)
#      my_rad_flx = read_my_array(profiles_file)
#  #    my_nc_up = read_my_array(profiles_file)
#  #    my_ss_up = read_my_array(profiles_file)
#  
#      print 'mean nc in cloud cells: ' , np.mean(my_nc[my_nc>20])
#    
#      profiles_file.close()
#      
#    
#      linestyles = ['--', '-.', ':']
#      dashList = [(3,1),(1,1),(4,1,1,1),(4,2)] 
#      if var == "thetal":
#        plot_my_array(axarr, plot_iter, my_thl, my_pos, nploty, xlabel=r'$\theta_l$ [K]', ylabel=ylabel, varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "qt":
#        plot_my_array(axarr, plot_iter, my_rtot, my_pos, nploty, xlabel='$q_{t}$ [g/kg]', ylabel=ylabel, varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "ql":
#        plot_my_array(axarr, plot_iter, my_rliq, my_pos, nploty, xlabel='$q_{l}$ [g/kg]', ylabel=ylabel, varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "cfrac":
#        plot_my_array(axarr, plot_iter, my_clfrac, my_pos, nploty, xlabel='Cloud fraction', ylabel=ylabel, varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "precip":
#        plot_my_array(axarr, plot_iter, my_prflux, my_pos, nploty, xlabel='Precip. flux [W m$^{-2}$]', ylabel=ylabel, varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "w_var":
#        plot_my_array(axarr, plot_iter, my_wvar, my_pos, nploty, xlabel=r'Var$\left(w\right)$ [m$^2$ s$^{-2}$]', ylabel=ylabel, varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "w_skw":
#        plot_my_array(axarr, plot_iter, my_w3rd, my_pos, nploty, xlabel='Third mom. of $w$ [m$^3$ s$^{-3}$]', ylabel=ylabel, varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "ss":
#        plot_my_array(axarr, plot_iter, my_ss, my_pos, nploty, xlabel='Supersaturation [\%]', ylabel=ylabel, varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "cl_nc":
#        plot_my_array(axarr, plot_iter, my_nc, my_pos, nploty, xlabel='$N_c$ [cm$^{-3}$]', ylabel=ylabel, varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "rad_flx":
#        plot_my_array(axarr, plot_iter, my_rad_flx, my_pos, nploty, xlabel='radiative flux [W m$^{-2}$]', ylabel=ylabel, varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "ndrop_cod_zoom":
#        plot_my_array(axarr, plot_iter, my_nc, my_pos, nploty, xlabel='$N_c$ [cm$^{-3}$]', ylabel=ylabel, varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#  #    # xrange of the nc_up plot
#  #    x = int(plot_iter / nploty)
#  #    y = plot_iter % nploty
#  #    axarr[x, y].set_ylim([0,1.2])
#  #    plot_iter = plot_my_array(axarr, plot_iter, my_nc_up, my_pos, xlabel='updraft nc', varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#  #    # xrange of the ss_up plot
#  #    x = int(plot_iter / nploty)
#  #    y = plot_iter % nploty
#  #    axarr[x, y].set_xlim([-5,1])
#  #    axarr[x, y].set_ylim([0,1.2])
#  #    plot_iter = plot_my_array(axarr, plot_iter, my_ss_up, my_pos, xlabel='updraft S', varlabel=profiles_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#  
#    except:
#      print 'error opening file: ', file_name
#      my_pos = 0
#      my_rtot = 0
#      my_rliq = 0
#      my_thl = 0
#      my_wvar = 0
#      my_w3rd = 0
#      my_prflux = 0
#      my_clfrac = 0
#      my_nc = 0
#      my_ss = 0
#      my_rad_flx = 0
#    label_counter = label_counter+1
#  plot_iter += 1
  dycoms_file.close()
  return plot_iter

def plot_series(var_list, plot_iter, nplotx, nploty, axarr, show_bin=False, suffix='', xlabel='', xlim=(1,6)):

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

  # create a list of files with UWLCM results
  series_files_names = []
  series_labels = []
  file_no = np.arange(1, len(sys.argv)-1 , 2)
  for no in file_no:
    series_files_names.append(argv[no] + suffix)
    series_labels.append(argv[no+1])

  for var in var_list:
    x = int(plot_iter / nploty)
    y = plot_iter % nploty

    # set plot limits
    if var == "surf_precip":
      axarr[x, y].set_ylim([-0.01,.5])

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
        
  #    axarr[x, y].fill_between(itime_h, minvar_arr, maxvar_arr, color='0.9')
  #    axarr[x, y].fill_between(itime_h, q1var_arr, q3var_arr, color='0.7')
  #    axarr[x, y].plot(itime_h, mvar_arr, color='black')
      axarr[x, y].set_xlim(xlim)
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
  
    
    
    
    #read my results

    
    label_counter=0
    for file_name in series_files_names:
      print file_name, var
      series_file = open(file_name, "r")

#      my_times = read_my_array(series_file)
#      my_max_w_var = read_my_array(series_file)
#      my_cfrac = read_my_array(series_file)
#      my_lwp = read_my_array(series_file)
#      my_er = read_my_array(series_file)
#      my_sp = read_my_array(series_file)
#      read_my_array(series_file) # discard accumulated precip
#      my_act_cond = read_my_array(series_file)
#      my_zb = read_my_array(series_file)
#      
#      # rescale time to hours
#      my_times = my_times / 3600.
#      
#      series_file.close()
#      
#  
#      linestyles = ['--', '-.', ':']
#      dashList = [(3,1),(1,1),(4,1,1,1),(4,2)] 
#      if var == "lwp":
#        plot_my_array(axarr, plot_iter, my_times, my_lwp, nploty, xlabel=xlabel, ylabel='LWP [g m$^{-2}$]', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "zi":
#        plot_my_array(axarr, plot_iter, my_times, my_er, nploty, xlabel=xlabel, ylabel='Entrainment rate [cm s$^{-1}$]', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "w2_max":
#        plot_my_array(axarr, plot_iter, my_times, my_max_w_var, nploty, xlabel=xlabel, ylabel='Max. $w$ variance [m$^{2}$ s$^{-2}$]', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "precip":
#        plot_my_array(axarr, plot_iter, my_times, my_sp, nploty, xlabel=xlabel, ylabel='Surface precip. [mm d$^{-1}$]', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "ndrop_cld":
#        plot_my_array(axarr, plot_iter, my_times, my_act_cond, nploty, xlabel=xlabel, ylabel='$N_c$ [cm$^{-3}$]', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#      if var == "zb":
#        plot_my_array(axarr, plot_iter, my_times, my_zb, nploty, xlabel=xlabel, ylabel='Cloud base height [m]', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])
#    # plot_my_array(axarr, plot_iter, my_times, my_cfrac, ylabel='Cloud fraction', varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])

      my_times = read_my_var(series_file, "position")
      my_res = read_my_var(series_file, var)

      # rescale time to hours
      my_times = my_times / 3600.

      series_file.close()

      linestyles = ['--', '-.', ':']
      dashList = [(3,1),(1,1),(4,1,1,1),(4,2)]
      plot_my_array(axarr, plot_iter, my_times, my_res, nploty, xlabel=xlabel, ylabel=var_labels[var], varlabel=series_labels[label_counter], dashes = dashList[label_counter % len(dashList)])

      label_counter+=1
  
    plot_iter = plot_iter + 1
  dycoms_file.close()
  return plot_iter
