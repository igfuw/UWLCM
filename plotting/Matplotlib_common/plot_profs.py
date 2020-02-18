import numpy as np
from sys import argv

from latex_labels import var_labels
from read_UWLCM_arrays import *

def plot_profiles(var_list, plot_iter, nplotx, nploty, axarr, xscaledict, yscaledict, xlimdict, ylimdict, suffix='', ylabel='', file_names=[], file_labels=[]):
  # if file names are not defined, read them and labels from command line
  if len(file_names)==0:
    file_no = np.arange(1, len(argv)-1 , 2)
    print file_no
    for no in file_no:
      file_names.append(argv[no]+suffix)
      file_labels.append(argv[no+1])
  
  for var in var_list:
    label_counter = 0
    for file_name in file_names:
      print file_name
    #  try:
      profiles_file = open(file_name, "r")
      my_pos = read_my_var(profiles_file, "position")
      my_res = read_my_var(profiles_file, var)
      profiles_file.close()

      # remove artificial values of cloudy-cell variables near the surface due to incorrect cloudiness mask (e.g. RICO mask used in DYCOMS)
      if(var == "gccn_rw_cl" or var == "non_gccn_rw_cl" or var == "cl_nc"):
        my_res[0:10] = 0
    
      linestyles = ['--', '-.', ':']
      dashList = [(3,1),(1,1),(4,1,1,1),(4,2)] 

      if(var == "base_prflux_vs_clhght" or var == "base_prflux_vs_clhght number of occurances"):
        plot_my_array(axarr, plot_iter, my_res, my_pos, nploty, xlabel=var_labels[var], ylabel="cloudy column height [m]", varlabel=file_labels[label_counter], dashes = dashList[label_counter % len(dashList)], xscale=xscaledict[var], yscale=yscaledict[var], xlim=xlimdict[var], ylim=ylimdict[var])
      else:
        plot_my_array(axarr, plot_iter, my_res, my_pos, nploty, xlabel=var_labels[var], ylabel=ylabel, varlabel=file_labels[label_counter], dashes = dashList[label_counter % len(dashList)], xscale=xscaledict[var], yscale=yscaledict[var], xlim=xlimdict[var], ylim=ylimdict[var])
   #   except:
   #     print 'error opening file: ', file_name
   #     my_pos = 0
   #     my_res = 0
      label_counter = label_counter+1
    plot_iter += 1
  return plot_iter

