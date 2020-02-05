from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../Matplotlib_common/")
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../Dycoms_RF02/")

from latex_labels import labeldict
from plot_ranges import * 
from plot_series import *
from plot_profs import *

# activate latex text rendering
rc('text', usetex=True)

dycoms_profs = ["prflux", "cl_nc", "non_gccn_rw_cl", "gccn_rw_cl"]
dycoms_series = ["surf_precip", "cl_gccn_conc"]
nplots = len(dycoms_profs + dycoms_series)# + 2 # 2 updraft profiles without dycoms results

# init the plot
nplotx = 2 #int(nplots/6 + 0.5)
nploty = 3 # int(float(nplots)/float(nplotx) + 0.5)
fig, axarr = plt.subplots(nplotx, nploty )

plot_iter=0
plot_iter = plot_series(dycoms_series, plot_iter, nplotx, nploty, axarr, xscaledict, yscaledict, xlimdict_series, ylimdict_series, False, suffix="series.dat", xlabel='Time [h]')
plot_profiles(dycoms_profs, plot_iter, nplotx, nploty, axarr, xscaledict, yscaledict, xlimdict_profs, ylimdict_profs, suffix="profiles_7200_18000.dat", ylabel='$z/z_i$')

# legend font size
plt.rcParams.update({'font.size': 8})

# hide axes on empty plots
if nplots % nploty == 0:
  nemptyplots=0
else:
  nemptyplots = nploty - nplots % nploty
emptyplots = np.arange(nploty - nemptyplots, nploty)
for empty in emptyplots:
  axarr[nplotx-1, empty].axis('off')

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
    # subplot numbering
    if y < nploty - nemptyplots or x < (nplotx - 1):
      axarr[x,y].text(0.8, 0.9, labeldict[y + x*nploty], fontsize=8, transform=axarr[x,y].transAxes)

## show legends
#for x in np.arange(nplotx):
#  for y in np.arange(nploty):
#    axarr[x,y].legend(loc="upper center")

#single legend for the whole figure
handles, labels = axarr[0,0].get_legend_handles_labels()
lgd = fig.legend(handles, labels, handlelength=4, loc='lower center', bbox_to_anchor=(0.45,0))


#figure size
fig.set_size_inches(7.874, 5. + (len(labels) - 2) * 0.2)# 5.214)#20.75,13.74)
#distances between subplots and from bottom of the plot
fig.subplots_adjust(bottom=0.14 + (len(labels) - 2) * 0.03, hspace=0.25, wspace=0.4)

#fig.tight_layout(pad=0, w_pad=0, h_pad=0)

#figure size
#fig.set_size_inches(7.874, 6 + (len(labels) - 2) * 0.2)# 5.214)#20.75,13.74)
#distances between subplots and from bottom of the plot
#fig.subplots_adjust(bottom=0.15 + (len(labels) - 2) * 0.02, hspace=0.25)


#plt.show()
fig.savefig(argv[len(sys.argv)-1], bbox_inches='tight', dpi=300)#, bbox_extra_artists=(lgd,))

