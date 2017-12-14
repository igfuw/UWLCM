#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

it_start = 0
it_end = 720
if len(sys.argv) == 3:
  it_start=sys.argv[1]
  it_end=sys.argv[2]

print it_start, it_end, len(sys.argv)

# load state
LoadState('paraview_DYCOMS_init_state.pvsm', LoadStateDataFileOptions='Use File Names From State',
    DataDirectory='/home/piotr/praca/wyniki/DYCOMS3D',
    OnlyUseFilesInDataDirectory=0,
    tempxmfFileNames=['/mnt/acc_prometheus/wyniki/DYCOMS3D/SD30_dt0.5_19_11_17/out_lgrngn/temp.xmf'],
    temp2xmfFileNames=['/mnt/acc_prometheus/wyniki/DYCOMS3D/SD30_dt1_4_12_17_outfreq30/out_lgrngn/temp.xmf'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2251, 1358]

# current camera placement for renderView1
#renderView1.CameraPosition = [3201.78803398526, -6348.11881882017, 4292.86270300142]
#renderView1.CameraFocalPoint = [3283.17952857712, 3129.45630056416, 514.117365215377]
#renderView1.CameraViewUp = [0.00, 0.370319799781735, 0.92889566640635]
#renderView1.CameraParallelScale = 4622.5

#renderView1.CameraPosition = [0, -0, 0]
#renderView1.CameraFocalPoint = [3200, 3200, 514]
#renderView1.CameraViewUp = [0., 0.8, 1]

#renderView1.CameraPosition = [3200, -6400, 750]
#renderView1.CameraFocalPoint = [3200, 320, 750]
#renderView1.CameraViewUp = [0.0, 0.0, 1.0]
#renderView1.CameraParallelScale = 4600

#camera = GetActiveCamera()
#camera.SetPosition(0,0,0) 
#camera.Elevation(8)



# find source
threshold6 = FindSource('Threshold6')
# set active source
SetActiveSource(threshold6)
# Adjust threshold and color map range for precipitation
threshold6.ThresholdRange = [0.005, 0.025]
# get color transfer function/color map for 'precip_rate'
precip_rateLUT = GetColorTransferFunction('precip_rate')
# Rescale transfer function
precip_rateLUT.RescaleTransferFunction(0.004, 0.025)
# get opacity transfer function/opacity map for 'precip_rate'
precip_ratePWF = GetOpacityTransferFunction('precip_rate')
# Rescale transfer function
precip_ratePWF.RescaleTransferFunction(0.004, 0.025)
# convert precip_rate color map from log to linear
precip_rateLUT.MapControlPointsToLinearSpace()
# Properties modified on precip_rateLUT
#precip_rateLUT.UseLogScale = 0
#display name of precip rate
# get color legend/bar for precip_rateLUT in view renderView1
precip_rateLUTColorBar = GetScalarBar(precip_rateLUT, renderView1)
# Properties modified on precip_rateLUTColorBar
precip_rateLUTColorBar.Title = 'precipitation flux [W / m^2]'
# invert coloring 
precip_rateLUT.InvertTransferFunction()


# hide color bar/color legend
threshold6Display = GetDisplayProperties(threshold6, view=renderView1)
threshold6Display.SetScalarBarVisibility(renderView1, False)
isovolume2 = FindSource('IsoVolume2')
isovolume2Display = GetDisplayProperties(isovolume2, view=renderView1)
isovolume2Display.SetScalarBarVisibility(renderView1, False)

# w color map range
SetActiveSource(isovolume2)
wLUT = GetColorTransferFunction('w')
wLUT.RescaleTransferFunction(-1.5, 2.)
# get opacity transfer function/opacity map for 'precip_rate'
wPWF = GetOpacityTransferFunction('w')
# Rescale transfer function
wPWF.RescaleTransferFunction(-1.5, 2.)


# add time display to the animation
# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter()
# Properties modified on annotateTimeFilter1
annotateTimeFilter1.Scale = 0.0002777777777778
# Properties modified on annotateTimeFilter1
annotateTimeFilter1.Format = 'Time: %.2f h'
# show data in view
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1)

# dont show domain borders
temp2xmf = FindSource('temp2.xmf')
# hide data in view
Hide(temp2xmf, renderView1)

# save animation
SaveAnimation('/home/piotr/praca/wyniki/DYCOMS3D/a02_SD30_dt1_4_12_17_outfreq30/paraview_batch_animation_'+it_start+'_'+it_end+'.avi', renderView1, ImageResolution=[560, 336],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0,
    ImageQuality=20,
    FrameRate=4,
    FrameWindow=[int(it_start), int(it_end)])

#### saving camera placements for all active views

# current camera placement for renderView1
#renderView1.CameraPosition = [3201.78803398526, -6348.11881882017, 4292.86270300142]
#renderView1.CameraFocalPoint = [3283.17952857712, 3129.45630056416, 514.117365215377]
#renderView1.CameraViewUp = [0.0040108379571942, 0.370319799781735, 0.92889566640635]
#renderView1.CameraParallelScale = 4622.5

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
