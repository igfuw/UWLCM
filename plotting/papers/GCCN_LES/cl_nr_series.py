from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../Matplotlib_common/")
#sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../Dycoms_RF02/")
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../RICO11/")

from plot_ranges import * 
from plot_series import *

# activate latex text rendering
rc('text', usetex=True)

series = ["cl_nr"]
nplotx = 1
nploty = 2

# init the plot
fig, axarr = plt.subplots(nplotx, nploty)
print 'axarr type: ', type(axarr)
print 'axarr shape: ', axarr.shape
plot_iter=0

assert len(argv) == 10


xlimdict_series = {
  0 : {"cl_nr" : (0,10)},
  1 : {"cl_nr" : (1,5)},
}


for cusc_iter in [0,1]: # cumulus, stratocumulus
  file_names = []
  file_labels = []
  file_no = np.arange(1 + 4 * cusc_iter, 1 + 4 * (cusc_iter+1) , 2)
  for no in file_no:
    print no
    print argv[no]
    file_names.append(argv[no] + "series.dat")
    file_labels.append(argv[no+1])
  
  print file_names
  plot_iter = plot_series(series, plot_iter, nplotx, nploty, axarr, xscaledict, yscaledict, xlimdict_series[cusc_iter], ylimdict_series, xlabel='Time [h]', file_names=file_names, file_labels=file_labels)

# legend font size
plt.rcParams.update({'font.size': 8})

#axes = plt.gca()
#axes.tick_params(direction='in')
labeldict=["{\it CuNc55}", "{\it ScNc45}"]
y_arr = np.arange(nploty)
for y in y_arr:
  #tics inside
  axarr[y].tick_params(direction='in', which='both', top=1, right=1)
  #minor tics
  axarr[y].xaxis.set_minor_locator(AutoMinorLocator())
  axarr[y].yaxis.set_minor_locator(AutoMinorLocator())
  #labels and tics font size
  for item in ([axarr[y].xaxis.label, axarr[y].yaxis.label] + axarr[y].get_xticklabels() + axarr[y].get_yticklabels()):
    item.set_fontsize(8)
  # subplot numbering
#    if y < nploty - nemptyplots or x < (nplotx - 1):
 #     axarr[y].text(0.8, 0.9, labeldict[y + x*nploty], fontsize=8, transform=axarr[y].transAxes)
  axarr[y].text(0.05, 0.85, labeldict[y], fontsize=8, transform=axarr[y].transAxes)

axarr[1].set_ylabel('')
axarr[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axarr[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axarr[1].xaxis.set_major_locator(MaxNLocator(integer=True))

## show legends
#for x in np.arange(nplotx):
#  for y in np.arange(nploty):
#    axarr[x,y].legend(loc="upper center")

#single legend for the whole figure
handles, labels = axarr[0].get_legend_handles_labels()
lgd = fig.legend(handles, labels, handlelength=4, loc='lower center', bbox_to_anchor=(0.45,0))


#figure size
fig.set_size_inches(7.874, 2. + (len(labels) ) * 0.34)# 5.214)#20.75,13.74)
#distances between subplots and from bottom of the plot
#fig.subplots_adjust(bottom=0. + (len(labels) ) * 0.044, hspace=0, wspace=0.4)

fig.tight_layout(pad=0, w_pad=1, h_pad=0, rect=(0,0.15,1,1))

#figure size
#fig.set_size_inches(7.874, 6 + (len(labels) - 2) * 0.2)# 5.214)#20.75,13.74)
#distances between subplots and from bottom of the plot
#fig.subplots_adjust(bottom=0.15 + (len(labels) - 2) * 0.02, hspace=0.25)

#hide vertical tick labels
#axarr[1].set_yticklabels([])
#axarr[2].set_yticklabels([])



#plt.show()
fig.savefig(argv[len(sys.argv)-1], bbox_inches='tight', dpi=300)#, bbox_extra_artists=(lgd,))

