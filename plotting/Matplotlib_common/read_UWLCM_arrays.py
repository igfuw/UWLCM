import numpy as np

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
  file_obj.seek(0)
  while True:
    arr, name = read_my_array(file_obj)
    if(str(name).strip() == str(var_name).strip()):
      break
  return arr

#def plot_my_array(axarr, plot_iter, time, val, nploty, xlabel=None, ylabel=None, varlabel=None , linestyle='--', dashes=(5,2), xlim=None, ylim=None, xscale="linear", yscale="linear"):
#  x = int(plot_iter / nploty)
#  y = plot_iter % nploty
#  if varlabel != None:
#    axarr[x, y].plot(time, val, label=varlabel, linestyle=linestyle, linewidth=1, dashes=dashes)
#  else:
#    axarr[x, y].plot(time, val, linestyle=linestyle, linewidth=1, dashes=dashes)
#  if xlabel:
#    axarr[x, y].set_xlabel(xlabel)
#  if ylabel:
#    axarr[x, y].set_ylabel(ylabel)
#  if xlim:
#    axarr[x, y].set_xlim(xlim)
#  if ylim:
#    axarr[x, y].set_ylim(ylim)
#  axarr[x, y].set_xscale(xscale)
#  axarr[x, y].set_yscale(yscale)

def plot_my_array(axarr, plot_iter, time, val, nploty, xlabel=None, ylabel=None, varlabel=None , linestyle='--', dashes=(5,2), xlim=None, ylim=None, xscale="linear", yscale="linear"):
  axflat = axarr.flatten() # flattening makes it work with 1d axarr
  if varlabel != None:
    axflat[plot_iter].plot(time, val, label=varlabel, linestyle=linestyle, linewidth=1, dashes=dashes)
  else:
    axflat[plot_iter].plot(time, val, linestyle=linestyle, linewidth=1, dashes=dashes)
  if xlabel:
    axflat[plot_iter].set_xlabel(xlabel)
  if ylabel:
    axflat[plot_iter].set_ylabel(ylabel)
  if xlim:
    axflat[plot_iter].set_xlim(xlim)
  if ylim:
    axflat[plot_iter].set_ylim(ylim)
  axflat[plot_iter].set_xscale(xscale)
  axflat[plot_iter].set_yscale(yscale)
