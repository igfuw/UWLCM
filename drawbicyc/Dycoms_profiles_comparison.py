#import h5py
from Dycoms_comparison_common import *

# activate latex text rendering
rc('text', usetex=True)

if sys.argv[len(sys.argv)-1] == "True":
  dycoms_vars = ["thetal", "qt", "ql", "cfrac", "precip", "w_var", "w_skw", "ss", "ndrop_cld", "rad_flx", "ndrop_cld_zoom"]
else:
  dycoms_vars = ["thetal", "qt", "ql", "cfrac", "precip", "w_var", "w_skw", "ss", "ndrop_cld", "rad_flx"]
nplots = len(dycoms_vars)# + 2 # 2 updraft profiles without dycoms results

# init the plot
nplotx = 2 #int(nplots/6 + 0.5)
nploty = int(float(nplots)/float(nplotx) + 0.5)
fig, axarr = plt.subplots(nplotx, nploty )

plot_iter=0
for var in dycoms_vars:
  if var == "thetal" or var == "w_var":
    plot_iter = plot_profiles(var, plot_iter, nplotx, nploty, axarr, ylabel='$z/z_i$')
  else:
    plot_iter = plot_profiles(var, plot_iter, nplotx, nploty, axarr)

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

# hide vertical tic labels
x_empty_label = np.arange(nplotx)
y_empty_label = np.arange(1, nploty)
for x in x_empty_label:
  for y in y_empty_label:
    axarr[x,y].set_yticklabels([])

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
fig.set_size_inches(7.874, 6 + (len(labels) - 2) * 0.2)# 5.214)#20.75,13.74)
#distances between subplots and from bottom of the plot
fig.subplots_adjust(bottom=0.18 + (len(labels) - 2) * 0.02, hspace=0.25)

#fig.tight_layout(pad=0, w_pad=0, h_pad=0)


#plt.show()
fig.savefig(argv[len(sys.argv)-2], bbox_inches='tight', dpi=300)#, bbox_extra_artists=(lgd,))

