from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../Matplotlib_common/")
from read_UWLCM_arrays import read_my_var

# activate latex text rendering
rc('text', usetex=True)

# assumed initial GCCN concentrations
GCCN_conc = [0,1,5,10]

# init the plot
nplotx = 2 #int(nplots/6 + 0.5)
fig, axarr = plt.subplots(1, nplotx )

#prepare a list of output files, assuming the following order: prsitine, standard, polluted, for each 4 results: no GCCN, GCCN ,GCCNx5, GCCNx10
file_names = []
file_no = np.arange(1, len(sys.argv)-1 , 1)
for no in file_no:
  file_names.append(sys.argv[no] + "series.dat")

assert(len(file_names) == 12)

label_counter=0

for it in np.arange(12):
  if(it % 4 == 0):
    mean_surf_precip = []
    tot_acc_surf_precip = []
  infile = open(file_names[it], "r")

# calc mean surf precip
  surf_precip = read_my_var(infile, "surf_precip")
  mean_surf_precip.append(np.mean(surf_precip[36:60])) # mean between 3rd and 5th hour. we assume outfreq of 300!
# calc acc surf precip
  acc_surf_precip = read_my_var(infile, "acc_precip")
  tot_acc_surf_precip.append(acc_surf_precip[60] - acc_surf_precip[12]) # acc precip between the 1st and the 5th hour, we assume outfreq of 300!
  if((it+1) % 4 == 0):
    #axarr[0].plot(GCCN_conc, mean_surf_precip, 'o')
    axarr[0].plot(GCCN_conc, tot_acc_surf_precip, 'o')



# legend font size
#plt.rcParams.update({'font.size': 8})

# hide axes on empty plots
#if nplots % nploty == 0:
#  nemptyplots=0
#else:
#  nemptyplots = nploty - nplots % nploty
#emptyplots = np.arange(nploty - nemptyplots, nploty)
#for empty in emptyplots:
#  axarr[nplotx-1, empty].axis('off')

#axes = plt.gca()
#axes.tick_params(direction='in')
#x_arr = np.arange(nplotx)
#y_arr = np.arange(nploty)
#for x in x_arr:
#  for y in y_arr:
#    #tics inside
#    axarr[x,y].tick_params(direction='in', which='both', top=1, right=1)
#    #minor tics
#    axarr[x,y].xaxis.set_minor_locator(AutoMinorLocator())
#    axarr[x,y].yaxis.set_minor_locator(AutoMinorLocator())
#    #labels and tics font size
#    for item in ([axarr[x,y].xaxis.label, axarr[x,y].yaxis.label] + axarr[x,y].get_xticklabels() + axarr[x,y].get_yticklabels()):
#      item.set_fontsize(8)
#    # subplot numbering
#    if y < nploty - nemptyplots or x < (nplotx - 1):
#      axarr[x,y].text(0.8, 0.9, labeldict[y + x*nploty], fontsize=8, transform=axarr[x,y].transAxes)

## show legends
#for x in np.arange(nplotx):
#  for y in np.arange(nploty):
#    axarr[x,y].legend(loc="upper center")

#single legend for the whole figure
#handles, labels = axarr[0,0].get_legend_handles_labels()
#lgd = fig.legend(handles, labels, handlelength=4, loc='lower center', bbox_to_anchor=(0.45,0))


#figure size
#fig.set_size_inches(7.874, 5. + (len(labels) - 2) * 0.2)# 5.214)#20.75,13.74)
#distances between subplots and from bottom of the plot
#fig.subplots_adjust(bottom=0.14 + (len(labels) - 2) * 0.03, hspace=0.25, wspace=0.4)

#fig.tight_layout(pad=0, w_pad=0, h_pad=0)

#figure size
#fig.set_size_inches(7.874, 6 + (len(labels) - 2) * 0.2)# 5.214)#20.75,13.74)
#distances between subplots and from bottom of the plot
#fig.subplots_adjust(bottom=0.15 + (len(labels) - 2) * 0.02, hspace=0.25)


plt.show()
#fig.savefig(argv[len(sys.argv)-1], bbox_inches='tight', dpi=300)#, bbox_extra_artists=(lgd,))

